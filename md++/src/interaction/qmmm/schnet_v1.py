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
#import schnetpack.transform as trn
#import pytorch_lightning as pl
import torch
#import torchmetrics
#import yaml

# ARGS                

# positions
# atomic_numbers
## device
## model_path
## val_model_paths
# step

# OUT

# energy
# forces
# val_dev

class SchNet_V1_Calculator:
    """
    A calculator class for molecular dynamics simulations using SchNetPack (spk) models. 

    This class provides tools to:
    - Predict energies and forces using a primary SchNet model.
    - Validate predictions against one or more additional validation models.
    - Compute the energy, forces, and validation deviations during simulations.

    Attributes:
        device (str): The device ('cuda' or 'cpu') to use for computations, set automatically.
        pred_calculator (spk.interfaces.SpkCalculator): The primary predictive calculator based on the provided model.
        val_calculators (list[spk.interfaces.SpkCalculator]): A list of validation calculators based on additional models.
        nn_valid_freq (int): Frequency (in time steps) at which validation predictions are performed.
        write_energy_freq (int): Frequency (in time steps) for writing energy predictions.

    Methods:
        get_environment(model_args):
            Retrieve the environment provider for a model based on its arguments.
        get_calculator(model_path):
            Initialize and return a SchnetPack calculator for a given model.
        predict_energy_and_forces(system):
            Predict the energy and forces for a given atomic system using the primary calculator.
        validate_prediction(system):
            Validate the energy predictions using additional validation calculators.
        calculate_next_step(atomic_numbers, positions, time_step):
            Compute energy, forces, and optionally validation deviation for a given atomic configuration.
    """

    def __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int) -> None:
        """
        Initialize the SchNet_V1_Calculator instance.

        Args:
            model_path (str): Path to the primary pre-trained model.
            val_model_paths (list): List of paths to additional validation models.
            nn_valid_freq (int): Frequency for validation predictions.
            write_energy_freq (int): Frequency for writing energy predictions.

        Sets up:
        - The primary predictive calculator (`pred_calculator`).
        - Optionally: Validation calculators (`val_calculators`) for models in `val_model_paths`.
        - Device (`device`) based on CUDA availability.
        """
        # Set computation device
        cuda_available = torch.cuda.is_available()
        self.device = 'cuda' if cuda_available else 'cpu'

        # Initialize the primary predictive calculator
        model_path = Path(model_path)
        self.pred_calculator = self.get_calculator(model_path=model_path)

        # Initialize validation calculators
        if len(val_model_paths) > 0:
            self.val_calculators = []
            for val_path in val_model_paths:
                val_path = Path(val_path)
                self.val_calculators.append(self.get_calculator(model_path=val_path))

        # Additional configuration
        self.nn_valid_freq = nn_valid_freq
        self.write_energy_freq = write_energy_freq
        return

    def get_environment(self, model_args) -> spk.environment.SimpleEnvironmentProvider:
        """
        Retrieve the environment provider for a machine learning model based on the specified arguments.

        Args:
            model_args (dict): A dictionary of model arguments required to configure the environment provider.

        Returns:
            spk.environment.SimpleEnvironmentProvider: An instance of the environment provider.
        """
        return spk.utils.script_utils.settings.get_environment_provider(model_args, device=self.device)

    def get_calculator(self, model_path) -> spk.interfaces.SpkCalculator:
        """
        Create and return a SchnetPack calculator instance using the specified model.

        Args:
            model_path (Path): Path to the pre-trained model file.

        Returns:
            spk.interfaces.SpkCalculator: A calculator configured with the model, environment provider, 
            and keys for energy and forces.
        """
        # Load the model and its configuration
        mlp_model = torch.load(model_path, map_location=self.device)
        model_args_path = model_path.parent / 'args.json'
        model_args = spk.utils.read_from_json(model_args_path.absolute())

        # Configure the environment provider
        mlp_environment = self.get_environment(model_args=model_args)

        # Create and return the calculator instance
        return spk.interfaces.SpkCalculator(
            model=mlp_model,
            device=self.device,
            environment_provider=mlp_environment,
            energy='energy',
            forces='forces'
        )

    def predict_energy_and_forces(self, system):
        """
        Predict the energy and forces for a given atomic system using the primary calculator.

        Args:
            system (ase.Atoms): An atomic system with positions and atomic numbers.

        Returns:
            tuple: Predicted energy (float) and forces (np.ndarray).
        """
        system.set_calculator(self.pred_calculator)
        energy = system.get_potential_energy()
        forces = system.get_forces()
        return energy, forces

    def validate_prediction(self, system):
        """
        Validate energy predictions using the additional validation calculators.

        Args:
            system (ase.Atoms): An atomic system with positions and atomic numbers.

        Returns:
            list: A list of energy predictions (float) from the validation calculators.
        """
        val_energies = []
        for val_calculator in self.val_calculators:
            system.set_calculator(val_calculator)
            val_energies.append(system.get_potential_energy())
        return val_energies

    def calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int)->None:
        """
        Compute the energy, forces, and optional validation deviation for a given atomic configuration.

        Args:
            atomic_numbers (list): Atomic numbers of the system's atoms.
            positions (list): Atomic positions in Cartesian coordinates.
            time_step (int): Current simulation time step.

        Returns:
            tuple: Energy (float), forces (np.ndarray), and optionally validation deviation (float).
        """
        # Create an ASE atomic system
        system = ase.Atoms(numbers=atomic_numbers, positions=positions)

        # Predict energy and forces
        self.energy, self.forces = self.predict_energy_and_forces(system=system)

        # Perform validation if applicable
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            val_energies = self.validate_prediction(system=system)
            if len(val_energies) == 1:
                self.nn_valid_dev = self.energy - val_energies[0]
            else:
                val_energies.append(self.energy)
                val_energies = np.array(val_energies)
                var = val_energies.var()
                self.nn_valid_dev = var / (len(val_energies) - 1)
        return None
    
    def get_energy(self):
        return self.energy
    
    def get_forces(self):
        return self.forces
    
    def get_nn_valid_dev(self):
        return self.nn_valid_dev

