"""
Example 6: Running Standard MD Simulation from Python
======================================================

This example shows how to run a standard molecular dynamics simulation
directly from Python using the GROMOS-RS MD engine.
"""

import gromos
import numpy as np
from pathlib import Path

print("=" * 70)
print("Example 6: Running Standard MD from Python")
print("=" * 70)

# NOTE: This example demonstrates the API.
# For actual simulations, you need:
# - Valid topology file (.top)
# - Valid coordinate file (.cnf)
# - Valid input parameter file (.imd)
# - Built GROMOS-RS binaries (cargo build --release)

print("""
To run MD simulations, you need the GROMOS-RS command-line binaries.

Build them with:
  cd gromos-rs
  cargo build --release

This will create:
  - gromos-rs/target/release/md       (standard MD)
  - gromos-rs/target/release/gamd     (GaMD)
  - gromos-rs/target/release/eds      (EDS)
  - gromos-rs/target/release/remd     (REMD)
""")

# ============================================================================
# Method 1: Using convenience function
# ============================================================================
print("\n" + "=" * 70)
print("Method 1: Quick Start with Convenience Function")
print("=" * 70)

print("""
from gromos import run_standard_md

# Run MD simulation
outputs = run_standard_md(
    topology="system.top",
    coordinates="start.cnf",
    input_file="md.imd",
    steps=10000,
    output_prefix="md"
)

print(f"Trajectory: {outputs['trajectory']}")
print(f"Energy: {outputs['energy']}")
print(f"Final config: {outputs['final_config']}")
""")

# ============================================================================
# Method 2: Using simulation class
# ============================================================================
print("\n" + "=" * 70)
print("Method 2: Using Simulation Class (More Control)")
print("=" * 70)

print("""
from gromos import MDSimulation

# Create simulation object
sim = MDSimulation(
    topology_file="system.top",
    coordinate_file="start.cnf",
    input_file="md.imd",
    output_prefix="production"
)

# Run simulation
outputs = sim.run(
    steps=100000,  # Override input file
    verbose=True   # Show output
)

# Access results
trajectory = outputs['trajectory']
energy = outputs['energy']
final_config = outputs['final_config']

# Analyze results (using NumPy/pandas/matplotlib)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read energy file (example - format depends on actual output)
# energies = pd.read_csv(energy, delim_whitespace=True)
# plt.plot(energies['time'], energies['total'])
# plt.xlabel('Time (ps)')
# plt.ylabel('Total Energy (kJ/mol)')
# plt.show()
""")

# ============================================================================
# Method 3: Creating a workflow
# ============================================================================
print("\n" + "=" * 70)
print("Method 3: Complete Workflow")
print("=" * 70)

print("""
def run_md_workflow(system_name, topology, coordinates, steps):
    '''
    Complete MD workflow:
    1. Energy minimization
    2. NVT equilibration
    3. NPT equilibration
    4. Production run
    '''
    from gromos import MDSimulation

    # 1. Energy minimization
    print("Step 1: Energy minimization...")
    em = MDSimulation(
        topology_file=topology,
        coordinate_file=coordinates,
        input_file="minimize.imd",
        output_prefix=f"{system_name}_em"
    )
    em_out = em.run(steps=1000)

    # 2. NVT equilibration
    print("Step 2: NVT equilibration...")
    nvt = MDSimulation(
        topology_file=topology,
        coordinate_file=str(em_out['final_config']),
        input_file="nvt.imd",
        output_prefix=f"{system_name}_nvt"
    )
    nvt_out = nvt.run(steps=10000)

    # 3. NPT equilibration
    print("Step 3: NPT equilibration...")
    npt = MDSimulation(
        topology_file=topology,
        coordinate_file=str(nvt_out['final_config']),
        input_file="npt.imd",
        output_prefix=f"{system_name}_npt"
    )
    npt_out = npt.run(steps=50000)

    # 4. Production run
    print("Step 4: Production MD...")
    prod = MDSimulation(
        topology_file=topology,
        coordinate_file=str(npt_out['final_config']),
        input_file="production.imd",
        output_prefix=f"{system_name}_prod"
    )
    prod_out = prod.run(steps=steps)

    return {
        'minimization': em_out,
        'nvt': nvt_out,
        'npt': npt_out,
        'production': prod_out
    }

# Run workflow
# results = run_md_workflow(
#     system_name="my_protein",
#     topology="protein.top",
#     coordinates="protein.cnf",
#     steps=1000000
# )
""")

# ============================================================================
# Python-based analysis
# ============================================================================
print("\n" + "=" * 70)
print("Post-Processing with Python")
print("=" * 70)

print("""
# After simulation, analyze results with Python tools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. Energy analysis
def analyze_energy(energy_file):
    '''Analyze energy trajectory'''
    # Read energy file
    data = pd.read_csv(energy_file, delim_whitespace=True)

    # Plot energies
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    axes[0, 0].plot(data['time'], data['kinetic'])
    axes[0, 0].set_title('Kinetic Energy')
    axes[0, 0].set_xlabel('Time (ps)')
    axes[0, 0].set_ylabel('Energy (kJ/mol)')

    axes[0, 1].plot(data['time'], data['potential'])
    axes[0, 1].set_title('Potential Energy')
    axes[0, 1].set_xlabel('Time (ps)')

    axes[1, 0].plot(data['time'], data['total'])
    axes[1, 0].set_title('Total Energy')
    axes[1, 0].set_xlabel('Time (ps)')
    axes[1, 0].set_ylabel('Energy (kJ/mol)')

    axes[1, 1].plot(data['time'], data['temperature'])
    axes[1, 1].set_title('Temperature')
    axes[1, 1].set_xlabel('Time (ps)')
    axes[1, 1].set_ylabel('Temperature (K)')

    plt.tight_layout()
    plt.savefig('energy_analysis.png', dpi=300)

    return data

# 2. Trajectory analysis
def analyze_trajectory(trajectory_file):
    '''Analyze trajectory'''
    # Read trajectory (format-specific)
    # ...

    # Calculate RMSD
    # Calculate RMSF
    # Calculate radius of gyration
    # etc.
    pass

# Run analysis
# energy_data = analyze_energy('md.tre')
# traj_data = analyze_trajectory('md.trc')
""")

# ============================================================================
# Integration with other tools
# ============================================================================
print("\n" + "=" * 70)
print("Integration with MDAnalysis / MDTraj")
print("=" * 70)

print("""
# GROMOS-RS outputs can be analyzed with standard tools

# Using MDAnalysis
import MDAnalysis as mda

# Load trajectory
u = mda.Universe('topology.pdb', 'trajectory.trc')

# Select atoms
protein = u.select_atoms('protein')
ligand = u.select_atoms('resname LIG')

# Calculate distances
from MDAnalysis.analysis import distances

# Calculate RMSD
from MDAnalysis.analysis.rms import RMSD
R = RMSD(protein, select='backbone')
R.run()

# Plot
import matplotlib.pyplot as plt
plt.plot(R.rmsd[:, 0], R.rmsd[:, 2])
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (Å)')
plt.savefig('rmsd.png')

# Using MDTraj
import mdtraj as md

# Load trajectory
traj = md.load('trajectory.trc', top='topology.pdb')

# Calculate RMSD
rmsd = md.rmsd(traj, traj, 0)

# Calculate contacts
distances, pairs = md.compute_contacts(traj)

# Secondary structure
dssp = md.compute_dssp(traj)
""")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("""
Running MD from Python:

1. Quick runs:
   from gromos import run_standard_md
   outputs = run_standard_md(topology, coords, input, steps=10000)

2. More control:
   from gromos import MDSimulation
   sim = MDSimulation(topology, coords, input)
   outputs = sim.run(steps=10000, verbose=True)

3. Complete workflows:
   - Chain multiple simulations
   - Minimization → equilibration → production
   - Use Python for orchestration

4. Analysis:
   - Use NumPy, pandas, matplotlib
   - Integrate with MDAnalysis, MDTraj
   - Full Python ecosystem available

Advantages:
✓ Pythonic interface
✓ Easy to script workflows
✓ Integration with analysis tools
✓ Rust performance under the hood
✓ Type safety and error handling

Next:
- See example 07 for GaMD
- See example 08 for EDS
- See example 09 for REMD
""")

print("\n" + "=" * 70)
print("Example complete!")
print("=" * 70)
