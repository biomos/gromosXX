import subprocess as sbp
import pandas as pd
import numpy as np
import os
from glob import glob
from importlib import reload
from shutil import rmtree
from pathlib import Path

import ase
import schnetpack as spk
import torch

class SchNet_V1_Calculator:
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
        energy (float): The predicted energy for the last system.
        forces (np.ndarray): The predicted forces for the last system.
        nn_valid_dev (float): Validation deviation calculated during validation.
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

    def get_environment(self, model_args) -> spk.environment.SimpleEnvironmentProvider:
        """
        Retrieve the environment provider for the model.

        Args:
            model_args (dict): Arguments to configure the environment provider.

        Returns:
            spk.environment.SimpleEnvironmentProvider: The configured environment provider.
        """
        return spk.utils.script_utils.settings.get_environment_provider(model_args, device=self.device)

    def get_calculator(self, model_path) -> spk.interfaces.SpkCalculator:
        """
        Create a SchnetPack calculator from the given model path.

        Args:
            model_path (Path): Path to the SchNet model.

        Returns:
            spk.interfaces.SpkCalculator: A configured SchnetPack calculator.
        """
        # Load the SchNet model and its configuration file.
        mlp_model = torch.load(model_path, map_location=self.device)
        model_args_path = model_path.parent / 'args.json'
        model_args = spk.utils.read_from_json(model_args_path.absolute())

        # Configure the environment provider.
        mlp_environment = self.get_environment(model_args=model_args)

        # Create and return the SchnetPack calculator.
        return spk.interfaces.SpkCalculator(
            model=mlp_model,
            device=self.device,
            environment_provider=mlp_environment,
            energy='energy',  # Energy key in model outputs.
            forces='forces'   # Forces key in model outputs.
        )

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
            if len(val_energies) == 1:
                self.nn_valid_dev = self.energy - val_energies[0]
            else:
                val_energies.append(self.energy)
                val_energies = np.array(val_energies)
                var = val_energies.var()
                self.nn_valid_dev = var / (len(val_energies) - 1)
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

    def get_nn_valid_dev(self):
        """
        Get the validation deviation calculated during the last validation step.

        Returns:
            float: Validation deviation.
        """
        return self.nn_valid_dev


