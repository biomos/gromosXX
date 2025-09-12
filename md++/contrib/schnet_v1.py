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

        calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int) -> None:
            Perform energy and force prediction, and validate the predictions if necessary.

        get_energy(self):
            Get the predicted energy for the last atomic system.

        get_forces(self):
            Get the predicted forces for the last atomic system.

        get_nn_valid_ene(self):
            Get the validation deviation calculated during the last validation step.
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

    def get_nn_valid_ene(self):
        """
        Get the validation deviation calculated during the last validation step.

        Returns:
            float: Validation deviation.
        """
        return self.nn_valid_dev

class Pert_SchNet_V1_Calculator_Error(Exception):
    pass

class Pert_SchNet_V1_Calculator(SchNet_V1_Calculator):
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
    def __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int, lam: float, perturbed_qm_states: list) -> None:
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
        super().__init__(model_path, val_model_paths, nn_valid_freq, write_energy_freq)
        
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
                # If state is B(2), add index to the 'B' list.
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
                dev = perturbed_energies_a_val.var() / (len(perturbed_energies_a_val) - 1)
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

        def calculate_forces(state, indices, lam_factor, shift: int = 0)->None:
            """
            Calculate the forces for a given state and set of indices.

            Args:
                state (int): The state for which to calculate the forces.
                indices (list of int): The indices of the particles for which to calculate the forces.
                lam_factor (float): The lambda factor used to interpolate between two sets of forces.
                shift (int, optional): The shift to apply to the indices when accessing the forces. Defaults to 0.

            Returns:
                None
            """
            for j, i in enumerate(indices):
                forces[i] = lam_factor * self.states[state]['forces_burnn'][j + shift] + (1 - lam_factor) * self.states[state]['forces_vac'][j + shift]
        
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
            for j, i in enumerate(indices):
                forces[i] = lam_factor * self.states['A']['forces_burnn'][j] + (1 - lam_factor) * self.states['A']['forces_vac'][j] + (1 - lam_factor) * self.states['B']['forces_burnn'][j] + lam_factor * self.states['B']['forces_vac'][j]

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
                calculate_forces('A', self.states_idx['only_A'], 1 - self.lam, shift=len(self.states_idx['both']))
                # Calculate forces for atoms only in state B.
                calculate_forces('B', self.states_idx['only_B'], self.lam, shift=len(self.states_idx['both']))
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

        # Calculate the perturbed energy and its derivative using the predicted energies.
        self.energy, self.derivative = self.calculate_perturbed_energy_and_derivative()
        self.forces = self.calculate_perturbed_forces()

        # Validate the perturbed energy if the current time step matches the validation frequency.
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            self.nn_valid_dev = self.validate_perturbed_energy()
        return None


    def get_derivative(self):
        """
        Get the derivative of the perturbed energy.

        Returns:
            float: Derivative of the perturbed energy.
        """
        return self.derivative
