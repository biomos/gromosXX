import numpy as np
from pathlib import Path
import os
import torch
import yaml
import ase
import schnetpack as spk

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
            energy_unit="eV", # units of energy property predicted by ase
            device= self.torchdevice, # device for computation
        )
        
        return calculator
        

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