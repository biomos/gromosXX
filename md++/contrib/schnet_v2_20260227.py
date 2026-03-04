import numpy as np
from pathlib import Path
import os
import torch
import yaml
import ase
import schnetpack as spk

EV_TO_KJMOL = 96.4853321233  
# Coulomb prefactor
K_E_KJMOL_NM_E2 = 138.935456  # kJ/mol·nm·e^-2
K_E_EV_NM_E2 = K_E_KJMOL_NM_E2 / EV_TO_KJMOL  # convert to eV·nm·e^-2

class YamlParser:
    def __init__(self, yaml_file):
        """
        Initializes the parser by loading the YAML file.
        
        :param yaml_file: Path to the YAML file.
        """
        self.data = self._load_yaml(yaml_file)
    
    def _load_yaml(self, yaml_file):
        """
        Loads and parses the YAML file.
        
        :param yaml_file: Path to the YAML file.
        :return: Parsed YAML data as a dictionary.
        """
        with open(yaml_file, 'r') as file:
            return yaml.safe_load(file)
    
    def get_cutoff(self):
        """
        Retrieves the cutoff value from the parsed YAML data.
        
        :return: The cutoff value from the YAML file.
        """
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
        """
        Attempts to identify the energy property key (usually starts with 'V' or contains 'ene').
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
        """
        Attempts to identify the force property key (usually starts with 'F' or contains 'force').
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
        """
        Identifies if the model was trained to predict charges.
        """
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
        """
        Identifies if the model was trained with electronic embedding.
        """
        try:
            representation = self.data['model']['representation']
            embeddings = representation.get('electronic_embeddings', [])
            if embeddings:
                return True
            return False
        except Exception:
            return False
 
class ExtendedConverter(spk.interfaces.AtomsConverter):
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

        # Add total_charge only if the model expects it
        if "total_charge" in atoms.info:
            charge = atoms.info.get("total_charge", 0.0)
            inputs["total_charge"] = torch.tensor([charge], dtype=torch.float32, device=self.device)

        # Add spin_multiplicity only if the model expects it
        if "spin_multiplicity" in atoms.info:
            multiplicity = atoms.info.get("spin_multiplicity", 0.0)
            inputs["spin_multiplicity"] = torch.tensor([multiplicity], dtype=torch.float32, device=self.device)
            
        # Add per-atom external potential phi only if present
        if "phi" in atoms.arrays:
            phi = atoms.arrays["phi"]
            phi = torch.tensor(phi, dtype=torch.float32, device=self.device)
            # ensure shape (N, 1) or (N,) depending on what your model expects
            inputs["phi"] = phi
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
        
        # Only provide charges_key if model was trained on partial charges
        charges_key = "charges" if self.has_partial_charges else None
        
        calculator = spk.interfaces.SpkCalculator(
            model_file = model_path, # path to model
            dtype=torch.float32, # 
            converter=ExtendedConverter,
            neighbor_list=spk.transform.ASENeighborList(cutoff=cutoff), # neighbor list
            energy_key=energy_key, # name of energy property in model
            force_key=force_key, # name of force property in model
            charges_key=charges_key,
            energy_unit="kJ/mol",
            position_unit="Angstrom",
            device= self.torchdevice, # device for computation
        )
        
        return calculator
        
    def _energy_forces_kjmol(self, system: ase.Atoms, calculator) -> tuple[float, np.ndarray]:
        system.calc = calculator
        energy_ev = system.get_potential_energy()      # ASE-native
        forces_ev_A = system.get_forces()              # ASE-native
        energy_kj = energy_ev * EV_TO_KJMOL
        forces_kj_A = forces_ev_A * EV_TO_KJMOL
        return energy_kj, forces_kj_A

    def predict_energy_and_forces(self, system: ase.Atoms):
        """
        Predict energy and forces for the given atomic system.

        Args:
            system (ase.Atoms): Atomic system for which predictions are required.

        Returns:
            tuple: Predicted energy (float) and forces (np.ndarray).
        """
        # Set the predictive calculator for the atomic system.
        system.calc = self.pred_calculator

        # Perform predictions.
        return self._energy_forces_kjmol(system=system, calculator=self.pred_calculator)
    
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
            e_kj, _ = self._energy_forces_kjmol(system=system, calculator=val_calculator)
            val_energies.append(e_kj)
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
            _, f_kj_A = self._energy_forces_kjmol(system=system, calculator=val_calculator)
            model_forces.append(np.linalg.norm(f_kj_A, axis=1))
            
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
        
        # per-atom potential
        #phi_kjmol = np.asarray(phi, dtype=np.float32)
        #phi_ev = phi_kjmol / EV_TO_KJMOL # internal ase units
        #system.arrays["phi"] = phi_ev

        # Predict energy and forces using the primary calculator.
        self.energy, self.forces = self.predict_energy_and_forces(system=system)
        
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

    def _coulomb_phi_ev(self, r_qeq_nm, r_or_nm, q_or_e, cutoff_nm):
        """
        Compute phi_i = k_e * sum_j q_j / r_ij (plain Coulomb with cutoff)
        Units:
            r in nm
            q in e
            returns phi in eV per e
        """

        if r_or_nm.shape[0] == 0:
            return torch.zeros(r_qeq_nm.shape[0], device=r_qeq_nm.device)

        # pair distances
        diff = r_qeq_nm[:, None, :] - r_or_nm[None, :, :]
        dist = torch.norm(diff, dim=-1)  # (Nq, Nor)

        # We dont need an extra cutoff we take the charge-group cutoff from gromos
        # mask = dist < cutoff_nm
        #safe_dist = torch.where(mask, dist, torch.ones_like(dist))
        #inv_r = torch.where(mask, 1.0 / safe_dist, torch.zeros_like(dist))
        inv_r = 1.0 / dist.clamp_min(1e-8)

        #phi = K_E_KJMOL_NM_E2 * torch.sum(inv_r * q_or_e[None, :], dim=1)

        phi = K_E_EV_NM_E2 * torch.sum(inv_r * q_or_e[None, :], dim=1)

        return phi

    def _model_forward_no_forces(self, inputs: dict) -> dict:
        """
        Forward through SchNetPack model while skipping any output module that
        computes forces via autograd. Still runs input_modules so PaiNN gets
        neighbor geometry (_Rij, idx_i/j, ...).
        """
        model = self.pred_calculator.model

        # 1) run input modules (neighbor list etc.)
        input_modules = getattr(model, "input_modules", None)
        if input_modules is not None:
            for mod in input_modules:
                inputs = mod(inputs)

        # 2) representation (PaiNN)
        inputs = model.representation(inputs)

        # 3) output modules (skip forces-producing ones)
        force_key = getattr(self.pred_calculator, "force_key", "forces")

        for mod in model.output_modules:
            outs = getattr(mod, "model_outputs", [])
            # Skip anything that outputs forces
            if (force_key in outs) or ("forces" in outs):
                continue
            inputs = mod(inputs)

        return inputs

    def calculate_next_step_b2(
        self,
        atomic_numbers: list,
        positions_A: list,
        time_step: int,
        or_positions_nm: list,
        or_charges_e: list,
        cutoff_nm: float,
    ) -> None:
        """
        B2: Compute phi internally from OR coords/charges, solve QEq in-model, and
        return conservative forces on both QEq sites and OR sites.
        Keeps same validation behavior as calculate_next_step().
        All returned forces are stored in kJ/mol/Å.
        """
        self.time_step = time_step

        # --- Build ASE system for validation (same as calculate_next_step) ---
        system = ase.Atoms(numbers=atomic_numbers, positions=positions_A)
        system.info["total_charge"] = float(self.total_charge)
        system.info["spin_multiplicity"] = float(self.spin_multiplicity)

        # --- Leaf coords (Å) for autograd ---
        r_qeq_A = torch.tensor(positions_A, dtype=torch.float32, device=self.torchdevice, requires_grad=True)

        if len(or_positions_nm) > 0:
            r_or_A = torch.tensor(or_positions_nm, dtype=torch.float32, device=self.torchdevice, requires_grad=True) * 10.0
            q_or = torch.tensor(or_charges_e, dtype=torch.float32, device=self.torchdevice)
        else:
            r_or_A = torch.zeros((0, 3), dtype=torch.float32, device=self.torchdevice, requires_grad=True)
            q_or = torch.zeros((0,), dtype=torch.float32, device=self.torchdevice)

        # nm versions only for Coulomb distance
        r_qeq_nm = r_qeq_A * 0.1
        r_or_nm  = r_or_A * 0.1

        # --- Prepare model inputs (converter) ---
        inputs = self.pred_calculator.converter(system)
        inputs[spk.properties.R] = r_qeq_A  # override with grad-enabled coords

        # --- Compute phi in Torch (eV per e) ---
        phi = self._coulomb_phi_ev(r_qeq_nm, r_or_nm, q_or, cutoff_nm)  # (Nq,)
        inputs["phi"] = phi

        # for validation models (they expect phi via ASE arrays)
        system.arrays["phi"] = phi.detach().cpu().numpy()

        # --- Forward model (must output charges + chi + Jii if using QEqCharges head) ---
        out = self._model_forward_no_forces(inputs)

        E_mlp = out[self.pred_calculator.energy_key].sum() / EV_TO_KJMOL # same unit as training

        # QEq outputs (your QEqCharges head) still not sure if units of chi and jii are correct they are not yet specified..the units are the same as 
        q   = out["charges"]
        chi = out["chi_qeq"]#*EV_TO_KJMOL
        Jii = out["Jii_qeq"]#*EV_TO_KJMOL
        
        # stationarity consistency check
        res = chi + phi + Jii * q  # all should be kJ/mol/e
        # per-molecule spread would be best, but global max spread is already informative
        print("stationarity spread (max abs dev from mean):",
            float((res - res.mean()).abs().max().detach().cpu()))
        print(chi, Jii,phi)
        # --- QEq energy (eV) ---
        E_qeq = torch.sum(chi * q + 0.5 * Jii * q**2 + q * phi)*E_mlp
        print(q,phi)
        E_total = E_mlp + E_qeq
        print(E_qeq,E_mlp,E_total)

        # --- Forces (eV/Å) ---
        dE_dRqeq, dE_dRor = torch.autograd.grad(
            E_mlp,
            [r_qeq_A, r_or_A],
            create_graph=False,
            retain_graph=False,
            allow_unused=True,
        )

        if dE_dRqeq is None:
            dE_dRqeq = torch.zeros_like(r_qeq_A)
        if dE_dRor is None:
            dE_dRor = torch.zeros_like(r_or_A)

        # Store (kJ/mol and kJ/mol/Å)
        self.energy = float(E_total.detach().cpu().item()) # if we fetch it with out[self.pred_calculator.energy_key].sum() we are already in kJ/mol
        self.forces = (-dE_dRqeq).detach().cpu().numpy()
        self.or_forces = (-dE_dRor).detach().cpu().numpy()

        # Optional: charges (same behavior)
        if self.has_partial_charges:
            self.charges = q.detach().cpu().numpy()

        # --- Validation (same logic as calculate_next_step) ---
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            val_energies = self.validate_prediction(system=system)
            self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(system=system)

            if len(val_energies) == 1:
                self.nn_valid_ene = (self.energy - val_energies[0]) / np.sqrt(2)
            else:
                val_energies = np.array(val_energies + [self.energy], dtype=float)
                self.nn_valid_ene = val_energies.std(ddof=1)

        return None

    def get_or_forces(self):
        return self.or_forces

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
    """
    Pert_SchNet_V1_Calculator is a specialized calculator for computing perturbed energies and forces 
    in molecular dynamics simulations using the SchNet model. It extends the SchNet_V1_Calculator class 
    to handle perturbations between different quantum mechanical (QM) end states.

    This class supports the following functionalities:
    - Initialization with model paths, validation frequencies, and perturbation settings.
    - Retrieval of atomic numbers and positions for specified QM states.
    - Calculation of perturbed energies and their derivatives.
    - Validation of perturbed energy calculations using multiple validation models.
    - Calculation of perturbed forces for the system.
    - Prediction of the next simulation step by updating energies, forces, and derivatives.

    Attributes:
        lam (float): Parameter for interpolating between different states.
        perturbed_qm_states (np.ndarray): Array specifying the QM states of each atom.
        qmzone_size (int): Length of the QM zone.
        num_states (int): Number of perturbed QM states.
        states_idx (dict): Dictionary storing indices of atoms in different states.
        states (dict): Dictionary storing atomic states and their properties.
        energy (float): Total perturbed energy.
        derivative (float): Derivative of the perturbed energy.
        forces (np.ndarray): Perturbed forces for the system.
        nn_valid_dev (float): Deviation of the perturbed energy during validation.

    Methods:
        __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int, lam: float, perturbed_qm_states: list) -> None:
            Initialize the Pert_SchNet_V1_Calculator with necessary model paths, validation frequency, and perturbation settings.

        get_z_and_positions_vac(self, atomic_numbers, positions, state_id):
            Retrieve the atomic numbers and positions within the quantum mechanical (QM) zone for a specified state.

        get_z_and_positions_vac_and_burnn(self, atomic_numbers, positions, state_id):
            Retrieve the atomic numbers and positions for both the vacuum and burnn configurations based on a specified state ID.

        get_state(self, atomic_numbers: list, positions: list, state_id: str) -> ase.Atoms:
            Construct two ASE `Atoms` objects representing the system end state in the vacuum and BuRNN environments for the specified perturbed quantum mechanical (QM) state.

        calculate_perturbed_energy_and_derivative(self) -> float:
            Calculate the perturbed energy and its derivative for two states (A and B) in both vacuum and Burnn environments.

        validate_state(self, state: str, calculator: spk.interfaces.SpkCalculator) -> float:
            Validate the state by calculating the potential energy using the specified validation calculator.

        validate_perturbed_energy(self):
            Validate the perturbed energy calculations using the provided end states in both vacuum and BuRNN environments.

        calculate_perturbed_forces(self):
            Calculate perturbed forces for the system.

        calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int) -> None:
            Calculate the next step in the MLP/MM simulation by updating the states, predicting energies and forces, and calculating the perturbed energy and forces.

        get_derivative(self):
            Get the derivative of the perturbed energy.
    """
    def __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int, 
                 spin_multiplicity:int,total_charge:int, lam: float, perturbed_qm_states: list) -> None:
        """
        Initialize the Pert_SchNet_V1_Calculator with necessary model paths, validation frequency, and perturbation settings.

        Args:
            model_path (str): Path to the primary SchNet model for energy and force calculations.
            val_model_paths (list): List of paths to validation models used to validate predictions.
            nn_valid_freq (int): Frequency (in steps) at which to validate the neural network predictions.
            write_energy_freq (int): Frequency (in steps) at which to write the energy predictions to output.
            lam (float): A parameter used for interpolating between different states.
            perturbed_qm_states (list): List specifying the QM states of each atom, where 0 indicates both states, 1 indicates state A, and 2 indicates state B.

        Returns:
            None
        """
        # Call the constructor of the parent class, SchNet_V1_Calculator, with the provided model paths and frequencies.
        super().__init__(model_path, val_model_paths, nn_valid_freq, write_energy_freq, spin_multiplicity, total_charge)
        
        # Store the perturbation parameter lambda (lam) for interpolating between states.
        self.lam = lam
        
        # Convert the list of perturbed QM states into a numpy array for easier manipulation and indexing.
        self.perturbed_qm_states = np.array(perturbed_qm_states)
        
        # Length of the QM zone.
        self.qmzone_size = len(self.perturbed_qm_states)

        # Initialize the `num_states` attribute with the number of perturbed QM states.
        self.num_states = len(set([s for s in self.perturbed_qm_states if s != 0]))

        # Create a dictionary to store indices of atoms in different states.
        self.states_idx = {}

        # Iterate over the perturbed QM states and populate the `states_idx` dictionary.
        for i, s in enumerate(self.perturbed_qm_states):
            if s == 0:
                # If state is both (0), add index to all relevant lists.
                try: 
                    self.states_idx['A'].append(i)
                    self.states_idx['B'].append(i)
                    self.states_idx['both'].append(i)
                except KeyError:
                    self.states_idx['A'] = [i]
                    self.states_idx['B'] = [i]
                    self.states_idx['both'] = [i]
            elif s == 1:
                # If state is 1, add index to the 'A' list.
                try: 
                    self.states_idx['A'].append(i)
                except KeyError: 
                    self.states_idx['A'] = [i]
            elif s == 2:
                # If state is 2, add index to the 'B' list.
                try: 
                    self.states_idx['B'].append(i)
                except KeyError: 
                    self.states_idx['B'] = [i]
            else:
                print("Error: Unknown state value:", s)
        
        # Check if both states are present and populate their respective lists.
        if 'both' in self.states_idx.keys():
            # Add indices of atoms only in A to the 'only_A' list.
            self.states_idx['only_A'] = [i for i in self.states_idx['A'] if i not in self.states_idx['both']]
            # Add indices of atoms only in B to the 'only_B' list.
            self.states_idx['only_B'] = [i for i in self.states_idx['B'] if i not in self.states_idx['both']]
        
        return


    def get_z_and_positions_vac(self, atomic_numbers, positions, state_id):
        """
        Retrieve the atomic numbers and positions within the quantum mechanical (QM) zone for a specified state.

        This method filters the atomic numbers and positions based on the provided state ID, returning
        only those that belong to the specified state within the QM zone.

        Args:
            atomic_numbers (list): A list of atomic numbers for the entire system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (str): The state identifier ('A' for state A, 'B' for state B) to filter the QM zone.

        Returns:
            tuple: Two lists containing:
                - atomic_numbers_vac (list): Atomic numbers within the QM zone for the specified state.
                - positions_vac (list): Positions corresponding to the atomic numbers within the QM zone.
        """
        atomic_numbers_vac = []  # Initialize an empty list to store atomic numbers for the specified state.
        positions_vac = []       # Initialize an empty list to store positions for the specified state.

        # Iterate over the indices of atoms in the specified state.
        for i in self.states_idx[state_id]:
            # Add atomic number and position to the respective lists.
            atomic_numbers_vac.append(atomic_numbers[i])
            positions_vac.append(positions[i])

        return atomic_numbers_vac, positions_vac


    def get_z_and_positions_vac_and_burnn(self, atomic_numbers, positions, state_id):
        """
        Retrieve the atomic numbers and positions for both the vacuum and burnn configurations based on a specified state ID.

        This method separates atomic numbers and positions into two groups:
        - Vacuum (vac): Atoms within the quantum mechanical (QM) zone corresponding to the specified state.
        - BuRNN (burnn): All atoms including the buffer zone.

        Args:
            atomic_numbers (list): A list of atomic numbers for the system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (str): The state identifier ('A' for state A, 'B' for state B) to filter the QM zone.

        Returns:
            tuple: Four lists containing:
                - atomic_numbers_vac (list): Atomic numbers within the QM zone for the specified state.
                - atomic_numbers_burnn (list): Atomic numbers for the entire Inner + Buffer zone, starting with those in the Inner zone.
                - positions_vac (list): Positions corresponding to the atomic numbers within the QM zone.
                - positions_burnn (list): Positions for the entire Inner + Buffer zone, starting with those in the Inner zone.
        """
        # Get atomic numbers and positions for the specified state within the QM zone.
        atomic_numbers_vac, positions_vac = self.get_z_and_positions_vac(
            atomic_numbers[:self.qmzone_size],  # Slice atomic numbers up to the QM zone size.
            positions[:self.qmzone_size],       # Slice positions up to the QM zone size.
            state_id                            # State identifier to filter the QM zone.
        )

        # Create a copy of the vacuum atomic numbers to start the burnn configuration.
        atomic_numbers_burnn = atomic_numbers_vac.copy()

        # Create a copy of the vacuum positions to start the burnn configuration.
        positions_burnn = positions_vac.copy()

        # Iterate over the atomic numbers and positions beyond the QM zone.
        for i, a in enumerate(atomic_numbers[self.qmzone_size:]):
            # Append the atomic number to the burnn configuration list.
            atomic_numbers_burnn.append(a)

            # Append the position to the burnn configuration list.
            positions_burnn.append(positions[i + self.qmzone_size])

        return atomic_numbers_vac, atomic_numbers_burnn, positions_vac, positions_burnn


    def get_state(self, atomic_numbers: list, positions: list, state_id: str) -> ase.Atoms:
        """
        Construct two ASE `Atoms` objects representing the system end state in the vacuum and BuRNN environments 
        for the specified perturbed quantum mechanical (QM) state.

        This method divides the atomic numbers and positions into two states: vacuum and Burnn, based on 
        the provided `state_id`, and then creates corresponding ASE `Atoms` objects for each state.

        Args:
            atomic_numbers (list): A list of atomic numbers for the entire system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (str): The state identifier ('A' for state A, 'B' for state B) used to filter and construct end states.

        Returns:
            tuple: A pair of ASE `Atoms` objects:
                - state_vac (ase.Atoms): The vacuum state representing the Z and positions for the given end state.
                - state_burnn (ase.Atoms): The BuRNN state representing the Z and positions for the given end state, including the buffer region.
        """
        # Retrieve atomic numbers and positions for both vacuum and Burnn states for the specified state ID.
        atomic_numbers_vac, atomic_numbers_burnn, positions_vac, positions_burnn = self.get_z_and_positions_vac_and_burnn(atomic_numbers, positions, state_id)

        # Create ASE Atoms object for the vacuum state.
        state_vac = ase.Atoms(numbers=atomic_numbers_vac, positions=positions_vac)

        # Create ASE Atoms object for the Burnn state.
        state_burnn = ase.Atoms(numbers=atomic_numbers_burnn, positions=positions_burnn)

        return state_vac, state_burnn


    def calculate_perturbed_energy_and_derivative(self) -> float:
        """
        Calculate the perturbed energy and its derivative for two states (A and B) in both vacuum and Burnn environments.

        This method computes a linear combination of the energies for state A and state B in both the vacuum 
        and Burnn conditions, using the lambda (`lam`) parameter. It also calculates the derivative of the 
        perturbed energy with respect to lambda.

        Returns:
            tuple: A tuple containing:
                - perturbed_energy (float): The total perturbed energy.
                - perturbed_energy_derivative (float): The derivative of the perturbed energy with respect to lambda.
        """
        # Iterate over each state to calculate the perturbed energy and its derivative.
        for state in self.states.keys():
            if state == 'A':
                # Calculate perturbed energy and derivative for state A.
                self.states[state]['perturbed_energy'] = (1 - self.lam) * self.states[state]['energy_burnn'] + self.lam * self.states[state]['energy_vac']
                self.states[state]['derivative'] = -self.states[state]['energy_burnn'] + self.states[state]['energy_vac']
            else:
                # Calculate perturbed energy and derivative for state B.
                self.states[state]['perturbed_energy'] = self.lam * self.states[state]['energy_burnn'] + (1 - self.lam) * self.states[state]['energy_vac']
                self.states[state]['derivative'] = self.states[state]['energy_burnn'] - self.states[state]['energy_vac']

        # Return the perturbed energies and derivatives for state A
        if len(self.states.keys()) == 1:
            return self.states['A']['perturbed_energy'], self.states['A']['derivative']
        # Sum the perturbed energies and derivatives of states A and B to get the total perturbed energy and derivative.
        else:
            return self.states['A']['perturbed_energy'] + self.states['B']['perturbed_energy'], self.states['A']['derivative'] + self.states['B']['derivative']

    def validate_state(self, state: str, calculator: spk.interfaces.SpkCalculator) -> float:
        """
        Validate the state by calculating the potential energy using the specified validation calculator.

        Args:
            state (str): The state to validate. Expected values are 'A' or other states.
            calculator (spk.interfaces.SpkCalculator): The calculator to set for the states.

        Returns:
            float: The calculated potential energy based on the state and lambda value.
        """
        self.states[state]['vac'].set_calculator(calculator)
        self.states[state]['burnn'].set_calculator(calculator)

        energy_burnn = self.states[state]['burnn'].get_potential_energy()
        energy_vac = self.states[state]['vac'].get_potential_energy()

        if state == 'A':
            return (1 - self.lam) * energy_burnn + self.lam * energy_vac
        else:
            return self.lam * energy_burnn + (1 - self.lam) * energy_vac
    
    def validate_perturbed_energy(self):
        """
        Validate the perturbed energy calculations using the provided end states in both vacuum and BuRNN environments.

        This method applies validation calculators to the provided atomic states to compute the perturbed energies 
        for both states A and B in vacuum and Burnn environments. It then calculates the deviation between the 
        previously computed perturbed energies and the validation energies, using either a single or multiple 
        validation calculators.

        Returns:
            float: The maximum absolute deviation between the computed perturbed energies and the validation energies.
        """
        # Case when only one state (A) is present.
        if len(self.states.keys()) == 1:
            # Case when only one validation calculator is provided.
            if len(self.val_calculators) == 1:
                # Compute validation perturbed energy for state A.
                perturbed_energy_a_val = self.validate_state('A', self.val_calculators[0])

                # Calculate deviation between computed and validation perturbed energy.
                dev = abs(self.states['A']['perturbed_energy'] - perturbed_energy_a_val)
        
            else:
                # Initialize lists to store validation perturbed energies for multiple calculators.
                perturbed_energies_a_val = [self.states['A']['perturbed_energy']]
                for val_calculator in self.val_calculators:
                    # Compute validation perturbed energy for state A.
                    perturbed_energy_a_val = self.validate_state('A', val_calculator)
                    # Append the calculated validation energy to the list.
                    perturbed_energies_a_val.append(perturbed_energy_a_val)
                
                # Convert the list of validation energies to a numpy array.
                perturbed_energies_a_val = np.array(perturbed_energies_a_val)
                # Calculate variance of the validation energies.
                dev = perturbed_energies_a_val.std()
            return dev

        # Case when two states (A and B) are present.
        else:
            # Case when only one validation calculator is provided.
            if len(self.val_calculators) == 1:
                # Compute validation perturbed energies for state A and state B.
                perturbed_energy_a_val = self.validate_state('A', self.val_calculators[0])
                perturbed_energy_b_val = self.validate_state('B', self.val_calculators[0])

                # Calculate deviations between computed and validation perturbed energies.
                dev_a = abs(self.states['A']['perturbed_energy'] - perturbed_energy_a_val)
                dev_b = abs(self.states['B']['perturbed_energy'] - perturbed_energy_b_val)

            else:
                # Initialize lists to store validation perturbed energies for multiple calculators.
                perturbed_energies_a_val = [self.states['A']['perturbed_energy']]
                perturbed_energies_b_val = [self.states['B']['perturbed_energy']]

                # Iterate over each validation calculator and compute perturbed energies.
                for val_calculator in self.val_calculators:
                    # Compute validation perturbed energies for state A and state B.
                    perturbed_energy_a_val = self.validate_state('A', val_calculator)
                    perturbed_energy_b_val = self.validate_state('B', val_calculator)

                    # Append the calculated validation energies to the lists.
                    perturbed_energies_a_val.append(perturbed_energy_a_val)
                    perturbed_energies_b_val.append(perturbed_energy_b_val)

                # Convert the lists of validation energies to numpy arrays.
                perturbed_energies_a_val = np.array(perturbed_energies_a_val)
                perturbed_energies_b_val = np.array(perturbed_energies_b_val)

                # Calculate variance of the validation energies for state A and state B.
                dev_a = perturbed_energies_a_val.var() / (len(perturbed_energies_a_val) - 1)
                dev_b = perturbed_energies_b_val.var() / (len(perturbed_energies_b_val) - 1)
            # Return the maximum absolute deviation between computed and validation energies.
            return np.max([dev_a, dev_b])


    def calculate_perturbed_forces(self):
        """
        Calculate the perturbed forces for the QM zone and buffer zone based on the lambda factor and the end states.
        This method calculates the perturbed forces for the QM zone and buffer zone by interpolating between the forces
        from different states ('A' and 'B') using the lambda factor. It handles cases where only one state ('A') is present
        or both states ('A' and 'B') are present.
        Returns:
            np.ndarray: A numpy array of shape (qm zone size + buffer size, 3) containing the calculated perturbed forces.
        """
        # Initialize the perturbed forces array, shape (qm zone size + buffer size, 3).
        forces = np.zeros((self.qmzone_size + (len(self.states['A']['forces_burnn']) - len(self.states['A']['forces_vac'])), 3), dtype=float)

        def calculate_forces(state, indices, lam_factor)->None:
            """
            Calculate the forces for a given state and set of indices.

            Args:
                state (int): The state for which to calculate the forces.
                indices (list of int): The indices of the particles for which to calculate the forces.
                lam_factor (float): The lambda factor used to interpolate between two sets of forces.

            Returns:
                None
            """
            for i in indices:
                i_state = self.states_idx[state].index(i)
                forces[i] = lam_factor * self.states[state]['forces_burnn'][i_state] + (1 - lam_factor) * self.states[state]['forces_vac'][i_state]
        
        def calculate_forces_both(indices, lam_factor)->None:
            """
            Calculate the combined forces for given indices using a lambda factor.

            This function computes the forces for each index in the provided list of indices
            by combining forces from different states ('A' and 'B') and conditions ('burnn' and 'vac')
            using the specified lambda factor.

            Args:
                indices (list of int): List of indices for which the forces need to be calculated.
                lam_factor (float): Lambda factor used to weight the contributions from different states and conditions.

            Returns:
                None: The function updates the `forces` array in place.
            """
            for i in indices:
                i_a = self.states_idx['A'].index(i)
                i_b = self.states_idx['B'].index(i)
                forces[i] = lam_factor * self.states['A']['forces_burnn'][i_a] + (1 - lam_factor) * self.states['A']['forces_vac'][i_a] + (1 - lam_factor) * self.states['B']['forces_burnn'][i_b] + lam_factor * self.states['B']['forces_vac'][i_b]

        # Case when only one state (A) is present.
        if len(self.states.keys()) == 1:
            # Calculate forces for the QM zone.
            forces[:self.qmzone_size] = (1 - self.lam) * self.states['A']['forces_burnn'][:self.qmzone_size] + self.lam * self.states['A']['forces_vac']
            # Calculate forces for the buffer zone.
            forces[self.qmzone_size:] = (1 - self.lam) * self.states['A']['forces_burnn'][self.qmzone_size:]
            return forces

        # Case when two states (A and B) are present.
        if len(self.states.keys()) == 2:
            if 'both' in self.states_idx.keys():
                # Calculate forces for atoms in both states.
                calculate_forces_both(self.states_idx['both'], 1 - self.lam)
                # Calculate forces for atoms only in state A.
                calculate_forces('A', self.states_idx['only_A'], 1 - self.lam)
                # Calculate forces for atoms only in state B.
                calculate_forces('B', self.states_idx['only_B'], self.lam)
            else:
                # Calculate forces for atoms in state A.
                calculate_forces('A', self.states_idx['A'], 1 - self.lam)
                # Calculate forces for atoms in state B.
                calculate_forces('B', self.states_idx['B'], self.lam)
            # Calculate forces for the buffer zone.
            forces[self.qmzone_size:] = (1 - self.lam) * self.states['A']['forces_burnn'][len(self.states_idx['A']):] + self.lam * self.states['B']['forces_burnn'][len(self.states_idx['B']):]
            return forces

    def calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int) -> None:
        """
        Calculate the next step in the MLP/MM simulation by updating the states, predicting energies and forces,
        and calculating the perturbed energy and forces.

        Args:
            atomic_numbers (list): List of atomic numbers for the atoms in the system.
            positions (list): List of positions for the atoms in the system.
            time_step (int): The current time step of the simulation.

        Returns:
            None
        """
        # Initialize the states dictionary to store vacuum and burnn states for each end state.
        self.states = {}

        # Iterate over the number of states (A and B) to populate the states dictionary.
        for i in range(self.num_states):
            if i == 0:
                state_id = 'A'
                self.states[state_id] = {'vac': [], 'burnn': []}
                # Retrieve the vacuum and burnn states for state A.
                self.states[state_id]['vac'], self.states[state_id]['burnn'] = self.get_state(atomic_numbers, positions, state_id=state_id)
            else:
                state_id = 'B'
                self.states[state_id] = {'vac': [], 'burnn': []}
                # Retrieve the vacuum and burnn states for state B.
                self.states[state_id]['vac'], self.states[state_id]['burnn'] = self.get_state(atomic_numbers, positions, state_id=state_id)

        # Predict energies and forces for each state in both vacuum and burnn environments.
        for state in self.states.keys():
            self.states[state]['energy_vac'], self.states[state]['forces_vac'] = self.predict_energy_and_forces(system=self.states[state]['vac'])
            self.states[state]['energy_burnn'], self.states[state]['forces_burnn'] = self.predict_energy_and_forces(self.states[state]['burnn'])
            print(f"{state}: forces Vac:\n{self.states[state]['forces_vac']}")
            print(f"{state}: forces BuRNN:\n{self.states[state]['forces_burnn']}")
        # Calculate the perturbed energy and its derivative using the predicted energies.
        self.energy, self.derivative = self.calculate_perturbed_energy_and_derivative()
        self.forces = self.calculate_perturbed_forces()
        print(self.forces)
        # Validate the perturbed energy if the current time step matches the validation frequency.
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            self.nn_valid_ene = self.validate_perturbed_energy()
        return None


    def get_derivative(self):
        """
        Get the derivative of the perturbed energy.

        Returns:
            float: Derivative of the perturbed energy.
        """
        return self.derivative
