#import subprocess as sbp
#import pandas as pd
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

class Pert_SchNet_V1_Calculator(SchNet_V1_Calculator):
    """
    A calculator class for perturbed energy and force calculations using modified SchNet models.

    This class extends the `SchNet_V1_Calculator` to handle perturbed quantum mechanical (QM) states,
    allowing the calculation of perturbed energies and forces in different states, and free derivative.

    Attributes:
        lam (float): A perturbation parameter for interpolating between states.
        perturbed_qm_states (np.ndarray): Array indicating QM states of atoms (1 for state A, 2 for state B).
        qmzone_size (int): Total number of QM states.
        lenght_state_a (int): Number of atoms in state A.
        lenght_state_b (int): Number of atoms in state B.
    """

    def __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int, lam: float, perturbed_qm_states: list) -> None:
        """
        Initialize the Pert_SchNet_V1_Calculator with necessary model paths, validation frequency, and perturbation settings.

        Args:
            model_path (str): Path to the primary SchNet model for energy and force calculations.
            val_model_paths (list): List of paths to validation models used to validate predictions.
            nn_valid_freq (int): Frequency (in steps) at which to validate the neural network predictions.
            write_energy_freq (int): Frequency (in steps) at which to write the energy predictions to output.
            lam (float): A perturbation parameter used for interpolating between different states.
            perturbed_qm_states (list): List specifying the QM states of each atom, where 1 indicates state A and 2 indicates state B.

        Returns:
            None
        """
        # Call the constructor of the parent class, SchNet_V1_Calculator, with the provided model paths and frequencies.
        super().__init__(model_path, val_model_paths, nn_valid_freq, write_energy_freq)
        
        # Store the perturbation parameter lambda (lam) for interpolating between states.
        self.lam = lam
        
        # Convert the list of perturbed QM states into a numpy array for easier manipulation and indexing.
        self.perturbed_qm_states = np.array(perturbed_qm_states)
        
        # Length of the QM zone.
        self.qmzone_size = len(self.perturbed_qm_states)
        
        # Determine and store the number of atoms in state A
        self.lenght_state_a = len(self.perturbed_qm_states[self.perturbed_qm_states == 1])
        
        # Determine and store the number of atoms in state B
        self.lenght_state_b = len(self.perturbed_qm_states[self.perturbed_qm_states == 2])
        
        return


    def get_z_and_positions_vac(self, atomic_numbers, positions, state_id):
        """
        Retrieve the atomic numbers and positions within the quantum mechanical (QM) zone for a specified end state.

        This method filters the atomic numbers and positions based on the provided end state ID, returning
        only those that belong to the specified end state within the QM zone.

        Args:
            atomic_numbers (list): A list of atomic numbers for the entire system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (int): The state identifier (1 for state A, 2 for state B) to filter the QM zone.

        Returns:
            tuple: Two lists containing:
                - atomic_numbers_vac (list): Atomic numbers within the QM zone for the specified state.
                - positions_vac (list): Positions corresponding to the atomic numbers within the QM zone.
        """
        atomic_numbers_vac = []  # Initialize an empty list to store atomic numbers for the specified state.
        positions_vac = []       # Initialize an empty list to store positions for the specified state.

        # Iterate over the atomic numbers and their indices.
        for i, a in enumerate(atomic_numbers):
            # Check if the current atom's end state matches the specified state ID.
            if self.perturbed_qm_states[i] == state_id:
                atomic_numbers_vac.append(a)  # Add atomic number if it matches the state.
                positions_vac.append(positions[i])  # Add position if it matches the state.

        return atomic_numbers_vac, positions_vac


    def get_z_and_positions_vac_and_burnn(self, atomic_numbers, positions, state_id):
        """
        Retrieve the atomic numbers and positions for both the vacuum and burnn configurations based on a specified end state ID.

        This method separates atomic numbers and positions into two groups:
        - Vacuum (vac): Atoms within the quantum mechanical (QM) zone corresponding to the specified state.
        - BuRNN (burnn): All atoms including the buffer zone.

        Args:
            atomic_numbers (list): A list of atomic numbers for the system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (int): The state identifier (1 for state A, 2 for state B) to filter the QM zone.

        Returns:
            Four lists containing:
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


    def get_state(self, atomic_numbers: list, positions: list, state_id: int) -> ase.Atoms:
        """
        Construct two ASE `Atoms` objects representing the system end state in the vacuum and BuRNN environments 
        for the specified quantum mechanical (QM) state.

        This method divides the atomic numbers and positions into two states: vacuum and Burnn, based on 
        the provided `state_id`, and then creates corresponding ASE `Atoms` objects for each state.

        Args:
            atomic_numbers (list): A list of atomic numbers for the entire system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (int): The state identifier (1 for state A, 2 for state B) used to filter and construct end states.

        Returns:
            tuple: A pair of ASE `Atoms` objects:
                - state_vac (ase.Atoms): The vacuum state representing the numbers and positions for the given end state.
                - state_burnn (ase.Atoms): The BuRNN state representing the numbers and positions for the given end state, including the buffer region.
        """
        # Retrieve atomic numbers and positions for both vacuum and Burnn states for the specified state ID.
        atomic_numbers_vac, atomic_numbers_burnn, positions_vac, positions_burnn = self.get_z_and_positions_vac_and_burnn(atomic_numbers, positions, state_id)

        state_vac = ase.Atoms(numbers=atomic_numbers_vac, positions=positions_vac)

        state_burnn = ase.Atoms(numbers=atomic_numbers_burnn, positions=positions_burnn)

        return state_vac, state_burnn


    def calculate_perturbed_energy_and_derivative(self, state_a_energy_vac: float, state_a_energy_burnn: float, state_b_energy_vac: float, state_b_energy_burnn: float) -> float:
        """
        Calculate the perturbed energy and its derivative for two states (A and B) in both vacuum and Burnn environments.

        This method computes a linear combination of the energies for state A and state B in both the vacuum 
        and Burnn conditions, using the lambda (`lam`) parameter. It also calculates the derivative of the 
        perturbed energy with respect to lambda.

        Args:
            state_a_energy_vac (float): Energy of state A in the vacuum
            state_a_energy_burnn (float): Energy of state A in the BuRNN scheme
            state_b_energy_vac (float): Energy of state B in the vacuum
            state_b_energy_burnn (float): Energy of state B in the BuRNN scheme

        Returns:
            tuple: A tuple containing:
                - perturbed_energy (float): The total perturbed energy.
                - perturbed_energy_derivative (float): The derivative of the perturbed energy with respect to lambda
        """
        # Calculate the perturbed energy for state A using a linear combination of vacuum and Burnn energies.
        self.perturbed_energy_a = (1 - self.lam) * state_a_energy_burnn + self.lam * state_a_energy_vac

        # Calculate the perturbed energy for state B using a linear combination of vacuum and Burnn energies.
        self.perturbed_energy_b = self.lam * state_b_energy_burnn + (1 - self.lam) * state_b_energy_vac

        # Sum the perturbed energies of states A and B to get the total perturbed energy.
        perturbed_energy = self.perturbed_energy_a + self.perturbed_energy_b

        # Calculate the derivative of the perturbed energy with respect to lambda
        perturbed_energy_derivative = -state_a_energy_burnn + state_a_energy_vac + state_b_energy_burnn - state_b_energy_vac

        return perturbed_energy, perturbed_energy_derivative


    def validate_perturbed_energy(self, state_a_vac: ase.Atoms, state_a_burnn: ase.Atoms, state_b_vac: ase.Atoms, state_b_burnn: ase.Atoms):
        """
        Validates the perturbed energy calculations using the provided end states in both vacuum and BuRNN environments.

        This method applies validation calculators to the provided atomic states to compute the perturbed energies 
        for both states A and B in vacuum and Burnn environments. It then calculates the deviation between the 
        previously computed perturbed energies and the validation energies, using either a single or multiple 
        validation calculators.

        Args:
            state_a_vac (ase.Atoms): The ASE Atoms object representing state A in the vacuum
            state_a_burnn (ase.Atoms): The ASE Atoms object representing state A in the BuRNN scheme
            state_b_vac (ase.Atoms): The ASE Atoms object representing state B in the vacuum
            state_b_burnn (ase.Atoms): The ASE Atoms object representing state B in the BuRNN scheme

        Returns:
            float: The maximum absolute deviation between the computed perturbed energies and the validation energies.
        """
        if len(self.val_calculators) == 1:
            # Set the single validation calculator for all states.
            state_a_vac.set_calculator(self.val_calculators[0])
            state_a_burnn.set_calculator(self.val_calculators[0])
            state_b_vac.set_calculator(self.val_calculators[0])
            state_b_burnn.set_calculator(self.val_calculators[0])

            # Compute validation perturbed energies for state A and state B.
            perturbed_energy_a_val = (1 - self.lam) * state_a_burnn.get_potential_energy() + self.lam * state_a_vac.get_potential_energy()
            perturbed_energy_b_val = self.lam * state_b_burnn.get_potential_energy() + (1 - self.lam) * state_b_vac.get_potential_energy()

            # Calculate deviations between computed and validation perturbed energies.
            dev_a = self.perturbed_energy_a - perturbed_energy_a_val
            dev_b = self.perturbed_energy_b - perturbed_energy_b_val

        else:
            # Initialize lists to store validation perturbed energies for multiple calculators.
            perturbed_energies_a_val = [self.perturbed_energy_a]
            perturbed_energies_b_val = [self.perturbed_energy_b]

            # Iterate over each validation calculator and compute perturbed energies.
            for val_calculator in self.val_calculators:
                state_a_vac.set_calculator(val_calculator)
                state_a_burnn.set_calculator(val_calculator)
                state_b_vac.set_calculator(val_calculator)
                state_b_burnn.set_calculator(val_calculator)

                perturbed_energy_a_val = (1 - self.lam) * state_a_burnn.get_potential_energy() + self.lam * state_a_vac.get_potential_energy()
                perturbed_energy_b_val = self.lam * state_b_burnn.get_potential_energy() + (1 - self.lam) * state_b_vac.get_potential_energy()

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
        return np.max([abs(dev_a), abs(dev_b)])


    def calculate_perturbed_forces(self, state_a_forces_vac, state_a_forces_burnn, state_b_forces_vac, state_b_forces_burnn):
        """
        Calculate perturbed forces for the system.

        Args:
            state_a_forces_vac (np.ndarray): Forces for state A in vacuum.
            state_a_forces_burnn (np.ndarray): Forces for state A in BuRNN scheme.
            state_b_forces_vac (np.ndarray): Forces for state B in vacuum.
            state_b_forces_burnn (np.ndarray): Forces for state B in BuRNN scheme.

        Returns:
            np.ndarray: Perturbed forces for the system.
        """

        # Initialize the perturbed forces array, shape (qm zone size, 3).
        forces = np.zeros((len(state_a_forces_vac) + len(state_b_forces_vac) + (len(state_a_forces_burnn) - len(state_a_forces_vac)), 3), dtype=float)

        # Calculate forces for state A.
        forces[:self.lenght_state_a] = (1 - self.lam) * state_a_forces_burnn[:self.lenght_state_a] + self.lam * state_a_forces_vac

        # Calculate forces for state B.
        forces[self.lenght_state_a:self.lenght_state_a + self.lenght_state_b] = self.lam * state_b_forces_burnn[:self.lenght_state_b] + (1 - self.lam) * state_b_forces_vac

        # Calculate forces for the buffer zone.
        forces[self.lenght_state_a + self.lenght_state_b:] = (1 - self.lam) * state_a_forces_burnn[self.lenght_state_a:] + self.lam * state_b_forces_burnn[self.lenght_state_b:]

        return forces

    def calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int) -> None:
        """
        Calculate the next step in the simulation by predicting energies and forces for perturbed states A and B, 
        and updating the perturbed energy, derivative, and forces.

        This function creates the atomic states for two perturbed states (A and B) in both vacuum and BuRNN environments. 
        It then predicts the energy and forces for these states and calculates the perturbed energy 
        and its derivative. If the time step matches the validation frequency, it validates the perturbed energy.

        Args:
            atomic_numbers (list): List of atomic numbers representing the atoms in the system.
            positions (list): List of atomic positions corresponding to the atomic numbers.
            time_step (int): The current time step of the MD simulation.

        Returns:
            None
        """
        # Retrieve states A in vacuum and BuRNN environments.
        state_a_vac, state_a_burnn = self.get_state(atomic_numbers, positions, state_id=1)
        
        # Predict energy and forces for state A in vacuum.
        state_a_energy_vac, state_a_forces_vac = self.predict_energy_and_forces(system=state_a_vac)
        
        # Predict energy and forces for state A in BuRNN scheme.
        state_a_energy_burnn, state_a_forces_burnn = self.predict_energy_and_forces(system=state_a_burnn)

        # Retrieve states B in vacuum and BuRNN environments.
        state_b_vac, state_b_burnn = self.get_state(atomic_numbers, positions, state_id=2)
        
        # Predict energy and forces for state B in vacuum.
        state_b_energy_vac, state_b_forces_vac = self.predict_energy_and_forces(system=state_b_vac)
        
        # Predict energy and forces for state B in BuRNN scheme.
        state_b_energy_burnn, state_b_forces_burnn = self.predict_energy_and_forces(system=state_b_burnn)

        # Calculate the perturbed energy and its derivative using the predicted energies.
        self.energy, self.derivative = self.calculate_perturbed_energy_and_derivative(state_a_energy_vac, state_a_energy_burnn, state_b_energy_vac, state_b_energy_burnn)
        
        # Calculate the perturbed forces using the predicted forces.
        self.forces = self.calculate_perturbed_forces(state_a_forces_vac, state_a_forces_burnn, state_b_forces_vac, state_b_forces_burnn)

        # Validate the perturbed energy if the current time step matches the validation frequency.
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            self.nn_valid_dev = self.validate_perturbed_energy(state_a_vac, state_a_burnn, state_b_vac, state_b_burnn)

        return None


    def get_derivative(self):
        """
        Get the derivative of the perturbed energy.

        Returns:
            float: Derivative of the perturbed energy.
        """
        return self.derivative
