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
        n_caps: int = 0,
    ) -> None:
        system = ase.Atoms(numbers=list(atomic_numbers), positions=list(positions))

        self.energy, self.forces = self.predict_energy_and_forces(system, self.pred_calculator)
        self.sextet_weight = 0.0

        if len(self.val_calculators) > 0 and int(time_step) % self.nn_valid_freq == 0:
            val_energies = self.validate_prediction(system)
            self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(system)
            if len(val_energies) == 1:
                self.nn_valid_ene = abs(self.energy - val_energies[0]) / np.sqrt(2.0)
            else:
                arr = np.array(val_energies + [self.energy], dtype=float)
                self.nn_valid_ene = float(arr.std(ddof=1))
        return None

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
