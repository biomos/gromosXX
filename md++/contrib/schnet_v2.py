from pathlib import Path
import os
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
import yaml
import ase
import schnetpack as spk

class YamlParser:
    """Small helper to read SchNetPack training configuration (``config.yaml``).

    The SchNetPack calculators need to know:
    - neighbor list cutoff
    - which property keys correspond to energy/forces
    - whether charges were trained/predicted
    - whether electronic embedding inputs are expected

    Parameters
    ----------
    yaml_file
        Path to the YAML configuration file written by the SchNetPack training run.
    """
    def __init__(self, yaml_file):
        """Initializes the parser by loading the YAML file."""
        self.data = self._load_yaml(yaml_file)
    
    def _load_yaml(self, yaml_file):
        """Load a YAML file into a dictionary."""
        with open(yaml_file, 'r') as file:
            return yaml.safe_load(file)
    
    def get_cutoff(self):
        """Return neighbor list cutoff (Angstrom) from the config file."""
        try:
            return self.data['globals']['cutoff']
        except KeyError as e:
            raise KeyError(f"Missing key in YAML file: {e}")
        
    def get_property_keys(self):
        """
        Returns all property keys defined under data.property_units.
        """
        try:
            return list(self.data["data"]["property_units"].keys())
        except KeyError as e:
            raise KeyError(f"Missing key in YAML file: {e}")

    def get_energy_key(self):
        """Heuristically identify the energy key.

        SchNetPack uses arbitrary property names; most training configs use
        something like ``V`` or ``energy``. We keep a conservative heuristic:

        - prefer keys starting with ``v`` (case-insensitive)
        - or keys containing ``ene`` (case-insensitive)
        """
        try:
            property_keys = self.get_property_keys()
            # Heuristic: look for a key that starts with 'V' or contains 'ene'
            for key in property_keys:
                if key.lower().startswith("v") or "ene" in key.lower():
                    return key
            raise ValueError("No energy key found in property_units.")
        except Exception as e:
            raise ValueError(f"Error finding energy key: {e}")

    def get_force_key(self):
        """Heuristically identify the forces key.

        - prefer keys starting with ``f`` (case-insensitive)
        - or keys containing ``force`` (case-insensitive)
        """
        try:
            property_keys = self.get_property_keys()
            for key in property_keys:
                if key.lower().startswith("f") or "force" in key.lower():
                    return key
            raise ValueError("No force key found in property_units.")
        except Exception as e:
            raise ValueError(f"Error finding force key: {e}")
        
    def get_charge_prediction(self):
        """Return True if the model task outputs include a charge prediction head."""
        try:
            outputs = self.data['task']['outputs']
            for out in outputs:
                name = out.get('name', '').lower()
                if 'charge' in name:
                    return True
            return False
        except Exception:
            return False

    def get_electronic_embedding(self):
        """Return True if the model representation declares electronic embeddings."""
        try:
            representation = self.data['model']['representation']
            embeddings = representation.get('electronic_embeddings', [])
            if embeddings:
                return True
            return False
        except Exception:
            return False
 
class ExtendedConverter(spk.interfaces.AtomsConverter):
    """
    AtomsConverter that optionally forwards global electronic state inputs.

    SchNetPack models may take additional (global) inputs such as
    ``total_charge`` and/or ``spin_multiplicity`` (e.g., for electronic embedding).

    This converter checks for these values in ``atoms.info`` and inserts them into
    the model inputs.

    Notes
    -----
    - Values are passed as tensors shaped ``(1,)`` on the same device as the model.
    - If a model does not use these keys, SchNetPack will ignore them safely.
      We still only add them when present in ``atoms.info`` to avoid surprising
      behavior for plain models.
    """
    def __init__(self, *args, **kwargs):
        """
        Initialize the converter.
        Detects which extra keys the model expects.
        """
        super().__init__(*args, **kwargs)
        # Detect which inputs the model actually expects
        self.model_input_keys = getattr(self, "model_input_keys", [])

    def __call__(self, atoms):
        inputs = super().__call__(atoms)
        
        # NOTE: we only add these keys if they exist in atoms.info. This allows
        # using the same model with/without electronic state conditioning.
        if "total_charge" in atoms.info:
            charge = atoms.info.get("total_charge", 0.0)
            inputs["total_charge"] = torch.tensor([charge], dtype=torch.float32, device=self.device)

        if "spin_multiplicity" in atoms.info:
            multiplicity = atoms.info.get("spin_multiplicity", 0.0)
            inputs["spin_multiplicity"] = torch.tensor([multiplicity], dtype=torch.float32, device=self.device)

        return inputs
    

class SchNet_V2_Calculator:
    """
    A calculator class for predicting energy and forces using SchNet models.

    This class provides functionality to load SchNet models, perform energy, 
    force and optionally partial charge predictions for atomic systems
    and validate predictions using additional validation models.

    Attributes:
        device (str): The computation device ('cuda' or 'cpu').
        pred_calculator (spk.interfaces.SpkCalculator): The primary predictive calculator.
        val_calculators (list): A list of validation calculators.
        nn_valid_freq (int): Frequency for performing validation (in terms of time steps).
        write_energy_freq (int): Frequency for writing energy values.
        energy (float): The predicted energy for the current system.
        forces (np.ndarray): The predicted forces for the current system.
        nn_valid_dev (float): Validation deviation calculated during validation.

    Methods:
        model_trained_on_partial_charges(model_path: str) -> bool:
            Checks if a SchNetPack model can predict partial charges
            
        __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int) -> None:
            Initialize the calculator with model paths and configuration.

        get_environment(self, model_args) -> spk.environment.SimpleEnvironmentProvider:
            Retrieve the environment provider for the model.

        get_calculator(self, model_path) -> spk.interfaces.SpkCalculator:
            Create a SchnetPack calculator from the given model path.

        predict_energy_and_forces(self, system: ase.Atoms):
            Predict energy and forces for the given atomic system.

        predict_partial_charges(self, system: ase.Atoms):
            Predict partial charges for the given atomic system.

        validate_prediction(self, system: ase.Atoms):
            Validate predictions using the validation calculators.
            
        validate_prediction_maxForceDeviation(self, system: ase.Atoms):
            Calculate Maximum force committee disagreement among all atoms in a structure.

        calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int) -> None:
            Perform energy and force prediction, and validate the predictions if necessary.

        get_energy(self):
            Get the predicted energy for the last atomic system.

        get_forces(self):
            Get the predicted forces for the last atomic system.

        get_nn_valid_ene(self):
            Get the validation deviation calculated during the last validation step.
            
        get_nn_valid_maxF(self):
            Get the maximum force committee disagreement among all atom during the last validation step.
    """
    
    @staticmethod
    def model_trained_on_partial_charges(model_path: str) -> bool:
        """
        Check if a SchNetPack model was trained on partial charges.
        """
        model_dir = Path(model_path).parent if Path(model_path).name == "best_model" else Path(model_path)
        config_path = model_dir / "config.yaml"
        
        return YamlParser(config_path).get_charge_prediction()
        
    @staticmethod
    def model_trained_with_electronic_embedding(model_path: str) -> bool:
        """
        Check if a SchNetPack model was trained with electronic embedding.
        """
        model_dir = Path(model_path).parent if Path(model_path).name == "best_model" else Path(model_path)
        config_path = model_dir / "config.yaml"
        return YamlParser(config_path).get_electronic_embedding()
        

    def __init__(self, model_path: str, val_model_paths: list,
                 nn_valid_freq: int, write_energy_freq: int,
                 spin_multiplicity:int,total_charge:int) -> None:
        """
        Initialize the calculator.
        """
        # Device
        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.torchdevice = torch.device(self.device)

        # Detect if the primary model was trained on charges
        self.has_partial_charges = self.model_trained_on_partial_charges(model_path)

        # Initialize primary calculator
        self.pred_calculator = self.get_calculator(model_path=Path(model_path))

        # Validation calculators
        self.val_calculators = []
        for val_path in val_model_paths:
            self.val_calculators.append(self.get_calculator(model_path=Path(val_path)))

        self.nn_valid_freq = nn_valid_freq
        self.write_energy_freq = write_energy_freq
        self.spin_multiplicity = spin_multiplicity
        self.total_charge = total_charge

    
    def get_calculator(self, model_path) -> spk.interfaces.SpkCalculator:
        """
        Create a SchnetPack calculator from the given model path.

        Args:
            model_path (Path): Path to the SchNet model.

        Returns:
            spk.interfaces.SpkCalculator: A configured SchnetPack calculator.
        """
        model_args_path = model_path.parent / 'config.yaml'
        model_path = os.path.join(model_path.parent,'best_model')

        # get cutoff and properties from configuration file
        yaml_parser = YamlParser(model_args_path)
        cutoff = yaml_parser.get_cutoff()
        energy_key = yaml_parser.get_energy_key()
        force_key = yaml_parser.get_force_key()
        
        # Decide per-model whether charges are available
        has_charges_local = yaml_parser.get_charge_prediction()
        charges_key = "charges" if has_charges_local else None
        
        calculator = spk.interfaces.SpkCalculator(
            model_file = model_path, # path to model
            dtype=torch.float32, # 
            converter=ExtendedConverter,
            neighbor_list=spk.transform.ASENeighborList(cutoff=cutoff), # neighbor list
            energy_key=energy_key, # name of energy property in model
            force_key=force_key, # name of force property in model
            charges_key=charges_key,
            energy_unit="eV", # units of energy property predicted by ase
            device= self.torchdevice, # device for computation
        )
        
        return calculator
        

    def predict_energy_and_forces(self, system: ase.Atoms, pred_calculator):
        """
        Predict energy and forces for the given atomic system.

        Args:
            system (ase.Atoms): Atomic system for which predictions are required.

        Returns:
            tuple: Predicted energy (float) and forces (np.ndarray).
        """
        # Set the predictive calculator for the atomic system.
        system.calc = pred_calculator
        # Perform predictions.
        energy = system.get_potential_energy()
        forces = system.get_forces()
        return energy, forces
    
    def predict_partial_charges(self, system: ase.Atoms):
        """
        Predict partial charges for the given atomic system.

        Args:
            system (ase.Atoms): Atomic system for which predictions are required.

        Returns:
            np.ndarray: Partial charges.
        """
        # Set the predictive calculator for the atomic system.
        system.calc = self.pred_calculator

        # Perform predictions.
        partial_charges = system.get_charges()
        return partial_charges
    
    def validate_prediction(self, system: ase.Atoms):
        """
        Validate predictions using the validation calculators.

        Args:
            system (ase.Atoms): Atomic system to validate.

        Returns:
            list: Validation energies computed by the validation calculators.
        """
        val_energies = []
        for val_calculator in self.val_calculators:
            system.calc = val_calculator
            val_energies.append(system.get_potential_energy())
        return val_energies
    
    def validate_prediction_maxForceDeviation(self, system: ase.Atoms):
        """
        Validate predictions using the validation calculators.
        Calculate Maximum force committee disagreement among all atoms in a structure 
        based on formula 1 in https://pubs.acs.org/doi/10.1021/acs.jctc.4c01382 
        Args:
            system (ase.Atoms): Atomic system to validate.

        Returns:
            list: Validation energies computed by the validation calculators.
        """
        model_forces = [np.linalg.norm(self.forces,axis=1)] # initialize production model
        # add validation model forces to model_forces
        for val_calculator in self.val_calculators:
            system.calc = val_calculator
            model_forces.append(np.linalg.norm(system.get_forces(),axis=1))
        # Calculate ensemble force 
        ensemble_Falpha = np.mean(model_forces, axis=0)  # Average force acting on every atom from all calculators

        # Calculate sigmaF_alpha
        Falpha_diff = [(np.asarray(force - ensemble_Falpha))**2 for force in model_forces]
        sigmaF_alpha = np.sqrt(np.mean(Falpha_diff, axis=0))

        # Store results
        maxForce_deviation = max(sigmaF_alpha)
        return maxForce_deviation

    def calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int) -> None:
        """
        Perform energy and force prediction, and validate the predictions if necessary.

        Args:
            atomic_numbers (list): List of atomic numbers for the atoms in the system.
            positions (list): List of atomic positions in the system.
            time_step (int): Current simulation time step.
        """
        # Create an ASE atomic system using atomic numbers and positions.
        system = ase.Atoms(numbers=atomic_numbers, positions=positions)

        system.info["total_charge"] = float(self.total_charge)
        system.info["spin_multiplicity"] = float(self.spin_multiplicity)
        
        # Predict energy and forces using the primary calculator.
        self.energy, self.forces = self.predict_energy_and_forces(system=system, pred_calculator=self.pred_calculator)
        
        # Predict partial charges if requested
        if self.has_partial_charges == True:
            self.charges = self.predict_partial_charges(system=system)

        # Perform validation at specified intervals.
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            val_energies = self.validate_prediction(system=system)

            # Calculate validation deviation based on validation results.
            self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(system=system)
            if len(val_energies) == 1:
                self.nn_valid_ene = (self.energy - val_energies[0])/np.sqrt(2) # devided by sqrt(2) to match the sample standard deviation bellow
            else:
                val_energies.append(self.energy)
                val_energies = np.array(val_energies)
                self.nn_valid_ene = val_energies.std(ddof=1) # to calculate sample standard deviation
        return None

    def get_energy(self):
        """
        Get the predicted energy for the last atomic system.

        Returns:
            float: Predicted energy.
        """
        return self.energy

    def get_forces(self):
        """
        Get the predicted forces for the last atomic system.

        Returns:
            np.ndarray: Predicted forces.
        """
        return self.forces

    def get_nn_valid_ene(self):
        """
        Get the validation energy deviation calculated during the last validation step.

        Returns:
            float: Validation energy deviation.
        """
        return self.nn_valid_ene

    def get_nn_valid_maxF(self):
        """
        Get the maximum force committee disagreement among all atom during the last validation step.

        Returns:
            float: Validation force deviation.
        """
        return self.nn_valid_maxF
     
    def get_charges(self):
        """
        Get the predicted partial charges for the last atomic system.

        Returns:
            np.ndarray: Predicted partial charges.
        """
        return self.charges

class Pert_SchNet_V2_Calculator(SchNet_V2_Calculator):
    """Thermodynamic integration / alchemical perturbation wrapper.

    This class evaluates two endstates (A and B) and interpolates them with λ.

    TI (linear endstate interpolation)
    ---------------------------------
    - Energy:
        ``E(λ) = (1-λ) E_A + λ E_B``
    - Forces:
        ``F(λ) = (1-λ) F_A + λ F_B``
    - Derivative:
        ``dE/dλ = E_B - E_A``

    Dual topology / state labels
    ----------------------------
    The per-atom state label array ``perturbed_qm_states`` is defined for the **IR block**
    (i.e., for the QM zone atoms only; BR and caps are handled separately):

    - ``0``: present/active in both endstates
    - ``1``: active only in A
    - ``2``: active only in B

    For each endstate we build two auxiliary systems:
    - **vac**: active IR subset only (no BR, no caps)
    - **burnn**: active IR subset + all BR + caps

    We then compute a corrected energy:

        ``E_corr = E_burnn - E_vac``

    If ``reference_energies`` are provided (default in this file), we use those instead of
    evaluating ``E_vac`` on the fly.

    Ordering contract
    -----------------
    Input ordering from C++ is:

        ``[ IR ][ BR ][ caps ]``

    burnn ordering (constructed here) is:

        ``[ active IR ][ all BR ][ caps ]``

    To return forces/charges back to the driver in the original ordering, we *scatter*
    burnn outputs back into the full-size arrays, leaving inactive IR atoms as zeros.
    """

    def __init__(
        self,
        model_path: str,
        val_model_paths: list,
        nn_valid_freq: int,
        write_energy_freq: int,
        spin_multiplicity: int,
        total_charge: int,
        lam: float,
        perturbed_qm_states: list,
        perturbed_total_charge: int,
        perturbed_spin_multiplicity: int,
        ref_vacA=None,
        ref_vacB=None,
    ) -> None:
        super().__init__(model_path, val_model_paths, nn_valid_freq, write_energy_freq,
                         spin_multiplicity, total_charge)

        self.lam = float(lam)
        
        self.ref_vacA = ref_vacA
        self.ref_vacB = ref_vacB
        # Reference energies for vacuum correction. If set to None, the code will
        # explicitly evaluate E_vac for each endstate and step.   
        if self.ref_vacA is not None and self.ref_vacB is not None:
            self.reference_energies = {"A": self.ref_vacA,"B": self.ref_vacB}
        else:
            self.reference_energies = None
        
        self.perturbed_qm_states = np.array(perturbed_qm_states, dtype=int)
        self.qmzone_size = len(self.perturbed_qm_states)
        
        self.nn_valid_ene = 0.0
        self.nn_valid_maxF = 0.0  # optional, currently not supported for perturbation code

        # Endstate-specific electronic settings
        self.charge_A = int(total_charge)
        self.spin_A = int(spin_multiplicity)
        self.charge_B = int(perturbed_total_charge)
        self.spin_B = int(perturbed_spin_multiplicity)

        # Active IR indices for each endstate (within IR block only)
        self.idx_A = [i for i, s in enumerate(self.perturbed_qm_states) if s in (0, 1)]
        self.idx_B = [i for i, s in enumerate(self.perturbed_qm_states) if s in (0, 2)]
        
        if len(self.idx_A) == 0 or len(self.idx_B) == 0:
            raise ValueError("Endstate A or B has zero active QM atoms; check perturbed_qm_states.")


        # Allow different prediction calculators per endstate (kept for future flexibility)
        self.pred_calculators = {"A": self.pred_calculator, "B": self.get_calculator(Path(model_path))}
        
        # Validation calculators for state A and B (one per validation model)
        self.val_calculators_A = self.val_calculators
        self.val_calculators_B = []

        for val_path in val_model_paths:
            self.val_calculators_B.append(self.get_calculator(Path(val_path)))

        self.energy = None
        self.forces = None
        self.derivative = None
        self.q_A: Optional[np.ndarray] = None
        self.q_B: Optional[np.ndarray] = None

    
    def _build_endstate_systems(
        self,
        atomic_numbers: Sequence[int],
        positions: Sequence[Sequence[float]],
        endstate: str,
        n_caps: int,
    ):
        """Construct the (vac, burnn) ASE systems for one endstate.

        Parameters
        ----------
        atomic_numbers
            Full list in ordering ``[IR][BR][caps]``.
        positions
            Full positions in ordering ``[IR][BR][caps]`` (Angstrom).
        endstate
            "A" or "B".
        n_caps
            Number of cap atoms at the end of the list.

        Returns
        -------
        vac : ase.Atoms
            Active IR-only system for vacuum reference.
        burnn : ase.Atoms
            Embedded system: (active IR + all BR) + caps.
        ir_idx : list[int]
            Global indices (in full ordering) of active IR atoms.
        br_idx : list[int]
            Global indices (in full ordering) of BR atoms.
        cap_start : int
            Index in full ordering where caps begin.
        """
        n_total = len(atomic_numbers)
        n_caps = int(n_caps)
        if n_caps < 0 or n_caps > n_total:
            raise RuntimeError(f"Invalid n_caps={n_caps} for n_total={n_total}")
        
        # IR size is fixed by perturbed_qm_states length
        n_ir = self.qmzone_size
        
        # caps are a contiguous tail
        cap_start = n_total - n_caps  # end of (IR+BR)
        if cap_start < n_ir:
            raise RuntimeError(f"cap_start={cap_start} < n_ir={n_ir} -> ordering mismatch")

        # Active IR indices and electronic settings per endstate
        if endstate == "A":
            ir_idx = [i for i, s in enumerate(self.perturbed_qm_states) if s in (0, 1)]
            charge, spin = self.charge_A, self.spin_A
        else:
            ir_idx = [i for i, s in enumerate(self.perturbed_qm_states) if s in (0, 2)]
            charge, spin = self.charge_B, self.spin_B

        # BR indices: always included in burnn, never in vac
        br_idx = list(range(n_ir, cap_start))
        qm_idx_burnn = ir_idx + br_idx

        # Build vac: active IR only
        z_vac = [atomic_numbers[i] for i in ir_idx]
        r_vac = [positions[i] for i in ir_idx]
        vac = ase.Atoms(numbers=z_vac, positions=r_vac)
        vac.info["total_charge"] = float(charge)
        vac.info["spin_multiplicity"] = float(spin)

        # Build burnn: (active IR + all BR) + caps
        z_qmbr = [atomic_numbers[i] for i in qm_idx_burnn]
        r_qmbr = [positions[i] for i in qm_idx_burnn]

        z_caps = list(atomic_numbers[cap_start:])
        r_caps = list(positions[cap_start:])

        burnn = ase.Atoms(numbers=z_qmbr + z_caps, positions=r_qmbr + r_caps)
        burnn.info["total_charge"] = float(charge)
        burnn.info["spin_multiplicity"] = float(spin)

        return vac, burnn, ir_idx, br_idx, cap_start

    def _scatter_burnn_vector_to_full(
        self,
        vec_burnn: np.ndarray,
        ir_idx: List[int],
        br_idx: List[int],
        cap_start: int,
        n_total: int,
    ) -> np.ndarray:
        """Scatter a burnn-ordered vector back to full ``[IR][BR][caps]`` ordering.

        burnn ordering:
            ``[active IR][all BR][caps]``

        full ordering:
            ``[IR][BR][caps]``

        Inactive IR atoms (not in ``ir_idx``) are left as zeros in the output.
        """
        
        vec_burnn = np.asarray(vec_burnn)

        n_ir_active = len(ir_idx)
        n_br = len(br_idx)
        n_caps = n_total - int(cap_start)

        expected = n_ir_active + n_br + n_caps
        if vec_burnn.shape[0] != expected:
            raise RuntimeError(
                f"burnn vector length mismatch: got {vec_burnn.shape[0]}, expected {expected} "
                f"(n_ir_active={n_ir_active}, n_br={n_br}, n_caps={n_caps})."
            )

        # Allocate output with same trailing dimensions
        out_shape = (n_total,) + vec_burnn.shape[1:]
        vec_full = np.zeros(out_shape, dtype=vec_burnn.dtype)

        # Slices in burnn
        v_ir = vec_burnn[:n_ir_active]
        v_br = vec_burnn[n_ir_active:n_ir_active + n_br]
        v_caps = vec_burnn[n_ir_active + n_br:]

        # Scatter active IR
        for local_i, global_i in enumerate(ir_idx):
            vec_full[global_i] = v_ir[local_i]

        # Scatter BR
        for local_i, global_i in enumerate(br_idx):
            vec_full[global_i] = v_br[local_i]

        # Scatter caps (caps are contiguous at the end in full indexing)
        vec_full[cap_start:] = v_caps

        return vec_full


    def _evaluate_corrected_endstate(
        self,
        atomic_numbers: Sequence[int],
        positions: Sequence[Sequence[float]],
        endstate: str,
        n_caps: int,
    ) -> Tuple[float, np.ndarray]:
        """Compute corrected (embedded - vacuum) energy and full-order forces.

        Notes
        -----
        The burnn prediction uses the compact ordering::

            [ active IR ][ all BR ][ caps ]

        while the full system ordering is::

            [ IR ][ BR ][ caps ]

        This method scatters the compact forces back into the full ordering.
        For consistency with charge handling, scattering is performed via
        :meth:`_scatter_burnn_vector_to_full` (applied per Cartesian component).
        """
        n_total = len(atomic_numbers)
        vac, burnn, ir_idx, br_idx, cap_start = self._build_endstate_systems(
            atomic_numbers, positions, endstate, n_caps
        )

        E_burnn, F_burnn = self.predict_energy_and_forces(burnn, pred_calculator=self.pred_calculators[endstate])

        if self.reference_energies is None:
            E_vac, F_vac = self.predict_energy_and_forces(vac, pred_calculator=self.pred_calculators[endstate])
        else:
            E_vac = float(self.reference_energies[endstate])
            F_vac = np.zeros((len(vac), 3), dtype=float)

        E_corr = E_burnn - E_vac

        # burnn forces: [active IR][all BR][caps]
        n_ir_active = len(ir_idx)
        n_br = len(br_idx)
        F_ir_active = F_burnn[:n_ir_active]
        F_br = F_burnn[n_ir_active:n_ir_active + n_br]
        F_caps = F_burnn[n_ir_active + n_br:]

        # subtract vacuum forces for IR (only defined on active IR subset)
        F_ir_corr = F_ir_active - F_vac

        # Build the compact "burnn corrected" force array with the same ordering
        # that _scatter_burnn_vector_to_full expects: [active IR][BR][caps].
        F_burnn_corr = np.concatenate([F_ir_corr, F_br, F_caps], axis=0)

        # Scatter each Cartesian component back into full ordering [IR][BR][caps].
        F_full_cols = [
            self._scatter_burnn_vector_to_full(
                F_burnn_corr[:, k], ir_idx, br_idx, cap_start, n_total
            )
            for k in range(3)
        ]
        F_full = np.stack(F_full_cols, axis=1)

        return float(E_corr), np.asarray(F_full)

    def _evaluate_endstate_charges(
        self,
        atomic_numbers: Sequence[int],
        positions: Sequence[Sequence[float]],
        endstate: str,
        n_caps: int,
    ) -> np.ndarray:
        """Predict burnn charges for one endstate and scatter back to full ordering."""
        n_total = len(atomic_numbers)
        _, burnn, ir_idx, br_idx, cap_start = self._build_endstate_systems(
            atomic_numbers, positions, endstate, n_caps
        )

        burnn.calc = self.pred_calculators[endstate]
        q_burnn = np.asarray(burnn.get_charges())

        q_full = self._scatter_burnn_vector_to_full(q_burnn, ir_idx, br_idx, cap_start, n_total)
        # Charges are 1D, keep dtype float
        return np.asarray(q_full, dtype=float).reshape((n_total,))

    
    def _evaluate_corrected_endstate_with_calculator(
        self,
        atomic_numbers: Sequence[int],
        positions: Sequence[Sequence[float]],
        endstate: str,
        val_calculator: spk.interfaces.SpkCalculator,
        n_caps: int,
    ) -> float:
        """Corrected endstate energy using a provided calculator (committee model)."""
        vac, burnn, _, _, _ = self._build_endstate_systems(atomic_numbers, positions, endstate, n_caps)

        E_burnn, _ = self.predict_energy_and_forces(burnn, pred_calculator=val_calculator)

        if self.reference_energies is None:
            E_vac, _ = self.predict_energy_and_forces(vac, pred_calculator=val_calculator)
        else:
            E_vac = float(self.reference_energies[endstate])

        return float(E_burnn - E_vac)

    def calculate_next_step(
        self,
        atomic_numbers: Sequence[int],
        positions: Sequence[Sequence[float]],
        time_step: int,
        n_caps: int = 0,
    ) -> None:
        """Run TI inference for a single step and update energy/forces/derivative.

        Parameters
        ----------
        atomic_numbers, positions
            Full system in ordering ``[IR][BR][caps]``.
        time_step
            MD step index (used for periodic committee validation).
        n_caps
            Number of cap atoms at the end of the list.

        Side effects
        ------------
        Updates:
        - ``self.energy`` and ``self.forces`` (interpolated)
        - ``self.derivative`` (dE/dλ)
        - ``self.charges`` (interpolated, if the model supports charge prediction)
        - ``self.nn_valid_ene`` (energy-only committee deviation, if enabled)
        """
        lam = float(self.lam)

        E_A, F_A = self._evaluate_corrected_endstate(atomic_numbers, positions, "A", n_caps)
        E_B, F_B = self._evaluate_corrected_endstate(atomic_numbers, positions, "B", n_caps)

        self.energy = (1.0 - lam) * E_A + lam * E_B
        self.forces = (1.0 - lam) * F_A + lam * F_B
        self.derivative = E_B - E_A
        # Optional Partial charges (if supported by the model)
        if getattr(self, "has_partial_charges", False):
            self.q_A = self._evaluate_endstate_charges(atomic_numbers, positions, "A",n_caps)
            self.q_B = self._evaluate_endstate_charges(atomic_numbers, positions, "B",n_caps)
            self.charges = (1.0 - lam) * self.q_A + lam * self.q_B
        
        # Energy-only committee validation
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            # production corrected endstates
            E_A, _ = self._evaluate_corrected_endstate(atomic_numbers, positions, "A",n_caps)
            E_B, _ = self._evaluate_corrected_endstate(atomic_numbers, positions, "B",n_caps)
            E_prod = (1.0 - lam) * E_A + lam * E_B
            
             # validation energies (same TI expression, but with each validation model)
            E_vals = []
            for val_calc_A,val_calc_B in zip(self.val_calculators_A,self.val_calculators_B):
                E_A_val = self._evaluate_corrected_endstate_with_calculator(atomic_numbers, positions, "A", val_calc_A,n_caps)
                E_B_val = self._evaluate_corrected_endstate_with_calculator(atomic_numbers, positions, "B", val_calc_B,n_caps)
                E_vals.append((1.0 - lam) * E_A_val + lam * E_B_val)

            # store committee deviation similar to base class
            if len(E_vals) == 1:
                self.nn_valid_ene = abs(E_prod - E_vals[0]) / np.sqrt(2)
            else:
                arr = np.array(E_vals + [E_prod], dtype=float)
                self.nn_valid_ene = arr.std(ddof=1)
        return None

    def get_derivative(self) -> float:
        """Return the last computed dE/dlambda."""
        return float(self.derivative)

    def get_charges_A(self) -> np.ndarray:
        """Return endstate-A charges (full ordering), if available."""
        if self.q_A is None:
            raise AttributeError("q_A not available (charges not predicted or not evaluated yet).")
        return np.asarray(self.q_A)

    def get_charges_B(self) -> np.ndarray:
        """Return endstate-B charges (full ordering), if available."""
        if self.q_B is None:
            raise AttributeError("q_B not available (charges not predicted or not evaluated yet).")
        return np.asarray(self.q_B)

