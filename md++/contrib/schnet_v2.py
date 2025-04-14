import numpy as np
from pathlib import Path
import os
import torch
import yaml
import ase
import schnetpack as spk

#import schnetpack.transform as trn
#from glob import glob
#from importlib import reload
#from shutil import rmtree

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
        Retrieves the cutoff value from the parsed YAML data.
        
        :return: The cutoff value from the YAML file.
        """
        try:
            return list(self.data['data']['property_units'].keys())
        except KeyError as e:
            raise KeyError(f"Missing key in YAML file: {e}")

class SchNet_V2_Calculator:
    """
    A calculator class for predicting energy and forces using SchNet models.

    This class provides functionality to load SchNet models, perform energy and 
    force predictions for atomic systems, and validate predictions using additional 
    validation models.

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
        __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int) -> None:
            Initialize the calculator with model paths and configuration.

        get_environment(self, model_args) -> spk.environment.SimpleEnvironmentProvider:
            Retrieve the environment provider for the model.

        get_calculator(self, model_path) -> spk.interfaces.SpkCalculator:
            Create a SchnetPack calculator from the given model path.

        predict_energy_and_forces(self, system: ase.Atoms):
            Predict energy and forces for the given atomic system.

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

    def __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int) -> None:
        """
        Initialize the calculator with model paths and configuration.

        Args:
            model_path (str): Path to the primary SchNet model.
            val_model_paths (list): List of paths to validation models.
            nn_valid_freq (int): Frequency to validate predictions.
            write_energy_freq (int): Frequency to write energy predictions.
        """
        # Set computation device to CUDA if available, otherwise fallback to CPU.
        cuda_available = torch.cuda.is_available()
        self.device = 'cuda' if cuda_available else 'cpu'
        self.torchdevice = torch.device(self.device)

        # Initialize the primary predictive calculator.
        model_path = Path(model_path)
        self.pred_calculator = self.get_calculator(model_path=model_path)

        # Initialize validation calculators if validation models are provided.
        self.val_calculators = []
        if len(val_model_paths) > 0:
            for val_path in val_model_paths:
                val_path = Path(val_path)
                self.val_calculators.append(self.get_calculator(model_path=val_path))

        # Store additional configuration parameters.
        self.nn_valid_freq = nn_valid_freq
        self.write_energy_freq = write_energy_freq
        return

    def get_calculator(self, model_path) -> spk.interfaces.SpkCalculator:
        """
        Create a SchnetPack calculator from the given model path.

        Args:
            model_path (Path): Path to the SchNet model.

        Returns:
            spk.interfaces.SpkCalculator: A configured SchnetPack calculator.
        """
        # get cutoff from configuration file
        model_args_path = model_path.parent / 'config.yaml'
        cutoff = YamlParser(model_args_path).get_cutoff()
        property_keys = YamlParser(model_args_path).get_property_keys()
        model_path = os.path.join(model_path.parent,'best_model')
        
        calculator = spk.interfaces.SpkCalculator(
            model_file = model_path, # path to model
            dtype=torch.float32, # 
            neighbor_list=spk.transform.ASENeighborList(cutoff=cutoff), # neighbor list
            energy_key=property_keys[0], # name of energy property in model
            force_key=property_keys[1], # name of force property in model
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
        system.set_calculator(self.pred_calculator)

        # Perform predictions.
        energy = system.get_potential_energy()
        forces = system.get_forces()
        return energy, forces

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
            system.set_calculator(val_calculator)
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
            system.set_calculator(val_calculator)
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

        # Predict energy and forces using the primary calculator.
        self.energy, self.forces = self.predict_energy_and_forces(system=system)

        # Perform validation at specified intervals.
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            val_energies = self.validate_prediction(system=system)

            # Calculate validation deviation based on validation results.
            self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(system=system)
            if len(val_energies) == 1:
                self.nn_valid_ene = self.energy - val_energies[0]
            else:
                val_energies.append(self.energy)
                val_energies = np.array(val_energies)
                self.nn_valid_ene = val_energies.std()
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