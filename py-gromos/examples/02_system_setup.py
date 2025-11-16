"""
Example 2: System Setup and Simulation Basics
==============================================

Demonstrates creating a molecular system and setting up simulation parameters.
"""

import gromos
import numpy as np

print("=" * 60)
print("GROMOS-RS: System Setup")
print("=" * 60)

# Create simulation box
print("\n--- Simulation Box ---")

# Vacuum (no periodic boundaries)
box_vacuum = gromos.Box.vacuum()
print(f"Vacuum box: {box_vacuum}")

# Rectangular periodic box (3x3x3 nm)
box_rect = gromos.Box.rectangular(3.0, 3.0, 3.0)
print(f"Rectangular box: {box_rect}")
print(f"Box dimensions: {box_rect.dimensions()}")
print(f"Box volume: {box_rect.volume():.2f} nm³")

# Larger box (5x5x5 nm)
box_large = gromos.Box.rectangular(5.0, 5.0, 5.0)
print(f"Large box: {box_large}")
print(f"Large box volume: {box_large.volume():.2f} nm³")

# Create system state
print("\n--- System State ---")
num_atoms = 1000
state = gromos.State(
    num_atoms=num_atoms,
    num_temp_groups=1,
    num_energy_groups=1
)
print(f"Created state with {state.num_atoms()} atoms")

# Initialize random positions in the box
print("\n--- Initialize Positions ---")
positions = np.random.uniform(0, 3.0, size=(num_atoms, 3)).astype(np.float32)
state.set_positions(positions)
print(f"Set positions: {positions.shape}")
print(f"First 3 positions:\n{positions[:3]}")

# Initialize velocities (Maxwell-Boltzmann at 300K)
print("\n--- Initialize Velocities ---")
# Simplified: should be mass-weighted for real simulation
temperature = 300.0  # K
kb = 0.00831446261815324  # kJ/(mol·K)
vel_scale = np.sqrt(kb * temperature)  # Simplified, needs mass
velocities = np.random.randn(num_atoms, 3).astype(np.float32) * vel_scale
state.set_velocities(velocities)
print(f"Set velocities at T={temperature}K")
print(f"Velocity scale: {vel_scale:.4f}")

# Retrieve data as NumPy arrays (zero-copy views)
print("\n--- Access Data ---")
pos_array = state.positions()
vel_array = state.velocities()
force_array = state.forces()

print(f"Positions shape: {pos_array.shape}")
print(f"Velocities shape: {vel_array.shape}")
print(f"Forces shape: {force_array.shape}")

# Energy tracking
print("\n--- Energy Tracking ---")
energy = gromos.Energy(num_temperature_groups=1, num_energy_groups=1)
print(f"{energy}")

# Manually set some energies (normally computed by force calculation)
energy.clear()
print(f"Cleared energy: total = {energy.total():.2f} kJ/mol")

print(f"Kinetic energy: {energy.kinetic:.2f} kJ/mol")
print(f"Potential energy: {energy.potential:.2f} kJ/mol")
print(f"Bond energy: {energy.bond:.2f} kJ/mol")
print(f"LJ energy: {energy.lj:.2f} kJ/mol")
print(f"Coulomb energy: {energy.coulomb:.2f} kJ/mol")

# Get energy as dictionary
energy_dict = energy.to_dict()
print(f"\nEnergy dictionary:")
for key, value in energy_dict.items():
    print(f"  {key}: {value:.2f} kJ/mol")

# Configuration (combines state + energy + topology)
print("\n--- Configuration ---")
config = gromos.Configuration(
    num_atoms=100,
    num_temp_groups=1,
    num_energy_groups=1
)
print(f"{config}")

current_state = config.current_state()
print(f"Current state: {current_state}")

current_energy = config.current_energy()
print(f"Current energy: {current_energy}")

# Topology
print("\n--- Topology ---")
topology = gromos.Topology()
print(f"{topology}")
print(f"Atoms: {topology.num_atoms()}")
print(f"Bonds: {topology.num_bonds()}")
print(f"Angles: {topology.num_angles()}")
print(f"Dihedrals: {topology.num_dihedrals()}")

print("\n" + "=" * 60)
print("System setup completed successfully!")
print("=" * 60)
