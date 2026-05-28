from pathlib import Path
from typing import Sequence

import numpy as np
import torch
import ase
from mace.calculators import MACECalculator


class MACE_GROMOS_Calculator:
    """Minimal MACE ASE wrapper for the GROMOS NN_Worker.

    Scope deliberately matches nn_model_type_standard only:
    - no dynamic partial charges
    - no perturbation / dE/dlambda
    - optional energy committee validation via val_model_paths
    - optional max-force committee disagreement, implemented like schnet_v2.py

    Units are expected to be ASE/MACE units on the Python side:
    positions in Angstrom, energy in eV, forces in eV/Angstrom.
    The C++ NN_Worker already performs the GROMOS <-> NN unit conversion.
    """

    def __init__(
        self,
        model_path: str,
        val_model_paths: list,
        nn_valid_freq: int,
        write_energy_freq: int,
        spin_multiplicity: int = 1,
        total_charge: int = 0,
    ) -> None:
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.default_dtype = "float32"

        self.pred_calculator = self.get_calculator(model_path)
        self.val_calculators = [self.get_calculator(path) for path in val_model_paths]

        self.nn_valid_freq = int(nn_valid_freq)
        self.write_energy_freq = int(write_energy_freq)
        self.spin_multiplicity = int(spin_multiplicity)
        self.total_charge = int(total_charge)

        self.energy = 0.0
        self.forces = np.zeros((0, 3), dtype=float)
        self.nn_valid_ene = 0.0
        self.nn_valid_maxF = 0.0
        self.sextet_weight = 0.0

        # Committee energies from the most recent standard BuRQM component calls.
        # These are kept in ASE/MACE energy units (eV).
        self._committee_energies_irbr: np.ndarray | None = None
        self._committee_energies_br: np.ndarray | None = None

    def get_calculator(self, model_path: str) -> MACECalculator:
        path = Path(model_path)
        # Accept either a direct model file or a training directory containing a common model filename.
        if path.is_dir():
            candidates = [
                path / "best_model.model",
                path / "model.model",
                path / "best_model.pt",
                path / "model.pt",
            ]
            for candidate in candidates:
                if candidate.exists():
                    path = candidate
                    break
            else:
                raise FileNotFoundError(
                    f"No MACE model file found in {path}. Expected one of: "
                    + ", ".join(c.name for c in candidates)
                )

        return MACECalculator(
            model_paths=str(path),
            device=self.device,
            #default_dtype=self.default_dtype,
        )

    @staticmethod
    def predict_energy_and_forces(system: ase.Atoms, calculator: MACECalculator):
        system.calc = calculator
        energy = system.get_potential_energy()
        forces = system.get_forces()
        return float(energy), np.asarray(forces, dtype=float)

    def validate_prediction(self, system: ase.Atoms):
        val_energies = []
        for val_calculator in self.val_calculators:
            e_val, _ = self.predict_energy_and_forces(system, val_calculator)
            val_energies.append(e_val)
        return val_energies

    @staticmethod
    def committee_std(committee_energies: Sequence[float]) -> float:
        arr = np.asarray(committee_energies, dtype=float)
        if arr.size <= 1:
            return 0.0
        if arr.size == 2:
            # Historical two-model convention: std([E_pred, E_val], ddof=1)
            # equals |E_pred - E_val| / sqrt(2).
            return float(abs(arr[0] - arr[1]) / np.sqrt(2.0))
        return float(arr.std(ddof=1))

    def validate_prediction_maxForceDeviation(self, system: ase.Atoms):
        model_forces = [np.linalg.norm(self.forces, axis=1)]
        for val_calculator in self.val_calculators:
            _, f_val = self.predict_energy_and_forces(system, val_calculator)
            model_forces.append(np.linalg.norm(f_val, axis=1))

        ensemble_falpha = np.mean(model_forces, axis=0)
        falpha_diff = [(np.asarray(force - ensemble_falpha)) ** 2 for force in model_forces]
        sigma_falpha = np.sqrt(np.mean(falpha_diff, axis=0))
        return float(max(sigma_falpha))

    def calculate_next_step(
        self,
        atomic_numbers: Sequence[int],
        positions: Sequence[Sequence[float]],
        time_step: int,
        total_charge: int | None = None,
        n_caps: int = 0,
        br_only: bool = False,
    ) -> None:
        if total_charge is not None:
            self.total_charge = int(total_charge)

        system = ase.Atoms(numbers=list(atomic_numbers), positions=list(positions))
        system.info["charge"] = self.total_charge

        self.energy, self.forces = self.predict_energy_and_forces(system, self.pred_calculator)
        self.sextet_weight = 0.0

        if len(self.val_calculators) > 0 and int(time_step) % self.nn_valid_freq == 0:
            val_energies = self.validate_prediction(system)
            self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(system)

            # Keep the prediction model first, followed by validation models.
            # This preserves model-to-model pairing between IR+BR and BR calls,
            # so the BuRQM uncertainty can be computed directly on
            # E_k(IR+BR) - E_k(BR).
            committee_energies = np.asarray([self.energy] + val_energies, dtype=float)
            self.nn_valid_ene = self.committee_std(committee_energies)

            if br_only:
                self._committee_energies_br = committee_energies
                self._append_nnvalidation_file("BR", int(time_step))
                if self._committee_energies_irbr is not None:
                    if self._committee_energies_irbr.shape == self._committee_energies_br.shape:
                        burqm_committee = self._committee_energies_irbr - self._committee_energies_br
                        self.nn_valid_ene = self.committee_std(burqm_committee)
                    else:
                        raise ValueError(
                            "IR+BR and BR validation committees have different sizes: "
                            f"{self._committee_energies_irbr.shape} vs {self._committee_energies_br.shape}"
                        )
            else:
                self._committee_energies_irbr = committee_energies
                self._append_nnvalidation_file("IRBR", int(time_step))
        return None
    
    def _append_nnvalidation_file(self, label: str, time_step: int) -> None:
        """Write per-evaluation validation so repeated NN calls do not overwrite each other.

        Values are written in Python/ASE units: energy in eV and max-force spread in eV/Angstrom.
        GROMOS still converts the scalar legacy getters to GROMOS units for the normal energy output.
        """
        path = Path(f"nnvalidation_{label}.dat")
        eV2kJpmol = 96.4869
        
        if not path.exists():
            path.write_text("# step nn_valid_ene_kJpmol nn_valid_maxF_kJpmol_per_Ang\n")
        with path.open("a") as handle:
            handle.write(f"{time_step:12d} {self.nn_valid_ene*eV2kJpmol:20.12e} {self.nn_valid_maxF*eV2kJpmol*10:20.12e}\n")

    def get_energy(self):
        return float(self.energy)

    def get_forces(self):
        return np.asarray(self.forces, dtype=float)

    def get_nn_valid_ene(self):
        return float(self.nn_valid_ene)

    def get_nn_valid_maxF(self):
        return float(self.nn_valid_maxF)

    def get_sextet_weight(self):
        return float(self.sextet_weight)
    
    def get_derivative(self):
        return float(0)
        