"""
Example 8: Running EDS (Enveloping Distribution Sampling) from Python
======================================================================

Enveloping Distribution Sampling uses a smoothed envelope potential to
enhance sampling of multiple end-states simultaneously.
"""

import gromos
import numpy as np

print("=" * 70)
print("Example 8: Running EDS from Python")
print("=" * 70)

# ============================================================================
# What is EDS?
# ============================================================================
print("\n" + "=" * 70)
print("What is EDS?")
print("=" * 70)

print("""
Enveloping Distribution Sampling (EDS) is a free energy method that:

1. Samples multiple end-states in a single simulation
2. Uses a smooth envelope potential
3. Calculates free energy differences
4. Better phase space overlap than standard FEP

Envelope potential:
  V_EDS = -s * ln(Σᵢ exp(-Vᵢ/s))

Where:
  - Vᵢ are the end-state potentials
  - s is the smoothness parameter (kJ/mol)
  - Larger s = smoother envelope = better sampling
  - Smaller s = closer to individual states

Key parameter:
  - s (smoothness): 0.5-2.0 kJ/mol (typical 1.0)
  - num_states: Number of end-states (2-10)

Applications:
  - Free energy calculations
  - Alchemical transformations
  - Ligand binding affinities
  - pKa predictions
  - Mutation effects
""")

# ============================================================================
# Method 1: Quick EDS run
# ============================================================================
print("\n" + "=" * 70)
print("Method 1: Quick EDS Simulation")
print("=" * 70)

print("""
from gromos import run_eds

# Run EDS simulation
outputs = run_eds(
    topology="system.top",
    coordinates="equilibrated.cnf",
    input_file="eds.imd",
    num_states=4,        # Number of end-states
    smoothness=1.0,      # Smoothness parameter (kJ/mol)
    steps=100000,
    output_prefix="eds"
)

print(f"Trajectory: {outputs['trajectory']}")
print(f"Energy: {outputs['energy']}")
print(f"Free energy: {outputs['free_energy']}")
print(f"State probabilities: {outputs['state_probabilities']}")
""")

# ============================================================================
# Method 2: Using EDS class
# ============================================================================
print("\n" + "=" * 70)
print("Method 2: EDS Class with Custom Parameters")
print("=" * 70)

print("""
from gromos import EDSSimulation

# Create EDS simulation
eds = EDSSimulation(
    topology_file="system.top",
    coordinate_file="start.cnf",
    input_file="eds.imd",
    output_prefix="eds_production",
    num_states=2,        # 2 end-states (state A and B)
    smoothness=1.0       # 1.0 kJ/mol smoothness
)

# Run EDS
outputs = eds.run(
    steps=500000,
    verbose=True
)

# Analyze results
free_energy_file = outputs['free_energy']
prob_file = outputs['state_probabilities']
""")

# ============================================================================
# Parameter selection
# ============================================================================
print("\n" + "=" * 70)
print("Selecting EDS Parameters")
print("=" * 70)

print("""
Choosing smoothness (s):
------------------------
- Too small (< 0.5 kJ/mol):
  • Poor sampling between states
  • Trajectories stuck in single state
  • Bad convergence

- Too large (> 2.0 kJ/mol):
  • Too much perturbation
  • Loss of information
  • Reduced accuracy

- Optimal (0.5-2.0 kJ/mol):
  • Good sampling of all states
  • Reasonable transitions
  • Balanced accuracy/sampling

Rule of thumb:
  s ≈ kT = 0.6 * T/300  kJ/mol
  At 300K: s ≈ 0.6 kJ/mol
  At 310K: s ≈ 0.62 kJ/mol

Practical values:
  - Ligand binding: s = 1.0-1.5 kJ/mol
  - Mutations: s = 0.5-1.0 kJ/mol
  - pKa calculations: s = 0.5 kJ/mol

Number of states:
-----------------
- 2 states: Simple A→B transformation
- 3-4 states: Intermediate states help convergence
- 5-10 states: Complex transformations
- More states = more sampling required

Simulation length:
------------------
Minimum: 100,000 steps
Typical: 500,000 - 5,000,000 steps
Rule: Run until ΔG converged to < 1 kJ/mol
""")

# ============================================================================
# Complete EDS workflow
# ============================================================================
print("\n" + "=" * 70)
print("Complete EDS Workflow")
print("=" * 70)

print("""
def run_eds_workflow(topology_A, topology_B, coordinates):
    '''
    Complete EDS workflow for free energy calculation:
    1. Equilibration
    2. EDS simulation
    3. Free energy analysis with WHAM/MBAR
    '''
    from gromos import MDSimulation, EDSSimulation
    import numpy as np

    # Step 1: Equilibrate each end-state
    print("Step 1: Equilibrating end-states...")

    # Equilibrate state A
    md_A = MDSimulation(
        topology_file=topology_A,
        coordinate_file=coordinates,
        input_file="equilibrate.imd",
        output_prefix="stateA_eq"
    )
    eq_A = md_A.run(steps=50000)

    # Equilibrate state B
    md_B = MDSimulation(
        topology_file=topology_B,
        coordinate_file=coordinates,
        input_file="equilibrate.imd",
        output_prefix="stateB_eq"
    )
    eq_B = md_B.run(steps=50000)

    # Step 2: Run EDS
    print("Step 2: Running EDS...")
    eds = EDSSimulation(
        topology_file=topology_A,  # Combined topology
        coordinate_file=str(eq_A['final_config']),
        input_file="eds.imd",
        output_prefix="eds",
        num_states=2,
        smoothness=1.0
    )
    eds_out = eds.run(steps=500000)

    # Step 3: Calculate free energy
    print("Step 3: Calculating free energy...")
    dG = calculate_free_energy(eds_out['free_energy'])

    print(f"Free energy difference: {dG:.2f} ± {dG_error:.2f} kJ/mol")

    return dG

def calculate_free_energy(dlg_file):
    '''Calculate free energy from EDS output'''
    import pandas as pd

    # Read free energy trajectory
    data = pd.read_csv(dlg_file, delim_whitespace=True)

    # Extract energies for each state
    V_R = data['V_stateA'].values  # Reference state
    V_i = [data[f'V_state{i}'].values for i in range(len(states))]

    # Calculate free energies using MBAR or WHAM
    # (simplified - actual implementation needs pymbar)
    dG = -kb * T * np.log(
        np.mean(np.exp(-V_i[1]/kb/T)) / np.mean(np.exp(-V_R/kb/T))
    )

    return dG

# Run workflow
# dG = run_eds_workflow('ligand_bound.top', 'ligand_free.top', 'start.cnf')
""")

# ============================================================================
# Free energy analysis
# ============================================================================
print("\n" + "=" * 70)
print("Free Energy Analysis with MBAR")
print("=" * 70)

print("""
EDS data can be analyzed with:
1. WHAM (Weighted Histogram Analysis Method)
2. MBAR (Multistate Bennett Acceptance Ratio)

Using pymbar for analysis:

def analyze_eds_with_mbar(dlg_file, num_states, temperature=300):
    '''Analyze EDS simulation with MBAR'''
    import numpy as np
    import pandas as pd
    from pymbar import MBAR

    # Constants
    kb = 0.00831446261815324  # kJ/(mol*K)
    beta = 1.0 / (kb * temperature)

    # Read EDS output
    data = pd.read_csv(dlg_file, delim_whitespace=True)

    # Extract energies for each state
    u_kn = np.zeros((num_states, len(data)))
    for i in range(num_states):
        u_kn[i, :] = beta * data[f'V_state{i}'].values

    # Number of samples from each state
    # (all from EDS ensemble)
    N_k = np.array([len(data)] + [0] * (num_states - 1))

    # Initialize MBAR
    mbar = MBAR(u_kn, N_k)

    # Calculate free energies
    results = mbar.getFreeEnergyDifferences()
    Delta_f_ij = results[0]  # Free energy differences
    dDelta_f_ij = results[1]  # Uncertainties

    # Print results
    print("Free Energy Differences (kJ/mol):")
    for i in range(num_states):
        for j in range(i+1, num_states):
            dG = Delta_f_ij[i, j] / beta
            dG_err = dDelta_f_ij[i, j] / beta
            print(f"  State {i} → State {j}: "
                  f"{dG:>8.2f} ± {dG_err:.2f} kJ/mol")

    return Delta_f_ij / beta, dDelta_f_ij / beta

# Analyze EDS results
# dG, dG_err = analyze_eds_with_mbar('eds.dlg', num_states=4, temperature=300)
""")

# ============================================================================
# EDS vs other free energy methods
# ============================================================================
print("\n" + "=" * 70)
print("EDS vs Other Free Energy Methods")
print("=" * 70)

print("""
Comparison:

Standard FEP (Free Energy Perturbation):
  Pros:
  + Simple theory
  + Well-established
  + Direct ΔG calculation

  Cons:
  - Poor phase space overlap
  - Many λ-windows needed
  - Slow convergence
  - Difficult for large perturbations

TI (Thermodynamic Integration):
  Pros:
  + Smooth free energy profile
  + Good for complex systems

  Cons:
  - Many λ-windows
  - Numerical integration errors
  - Slow convergence

EDS (Enveloping Distribution Sampling):
  Pros:
  + Single simulation for all states
  + Better phase space overlap
  + Faster convergence
  + Can handle large perturbations
  + Multiple states simultaneously

  Cons:
  - Parameter tuning (s)
  - More complex analysis
  - Requires MBAR/WHAM

When to use EDS:
  ✓ Large perturbations
  ✓ Multiple states to compare
  ✓ Poor overlap in standard FEP
  ✓ Limited computational resources
  ✗ Very small perturbations
  ✗ Simple two-state systems
""")

# ============================================================================
# Common applications
# ============================================================================
print("\n" + "=" * 70)
print("Common EDS Applications")
print("=" * 70)

print("""
1. Ligand Binding Free Energy:
   States: Ligand bound / Ligand unbound
   ΔG_bind = G_complex - G_protein - G_ligand

   Example:
     State A: Protein-ligand complex
     State B: Protein + ligand in solution

2. Relative Binding Free Energy:
   Compare two ligands binding to same protein
   ΔΔG = ΔG_bind(B) - ΔG_bind(A)

   Example:
     State A: Protein-ligandA
     State B: Protein-ligandB

3. Mutation Free Energy:
   Compare wild-type vs mutant protein
   ΔΔG = G_mutant - G_wildtype

   Example:
     State A: Wild-type protein
     State B: F123A mutant

4. pKa Calculations:
   Protonation state free energies
   pKa = pH + ΔG/(2.303 * RT)

   Example:
     State A: Protonated residue
     State B: Deprotonated residue

5. Conformational Free Energy:
   Different conformations of same molecule

   Example:
     State A: α-helix
     State B: β-sheet
""")

# ============================================================================
# Best practices
# ============================================================================
print("\n" + "=" * 70)
print("EDS Best Practices")
print("=" * 70)

print("""
1. System preparation:
   - Equilibrate each end-state separately
   - Ensure proper solvation
   - Check for steric clashes
   - Verify charge neutrality

2. Parameter selection:
   - Start with s = 1.0 kJ/mol
   - Adjust based on convergence
   - Test different values

3. Simulation length:
   - Minimum 100,000 steps
   - Monitor ΔG convergence
   - Run until error < 1 kJ/mol

4. Analysis:
   - Use MBAR (preferred) or WHAM
   - Calculate error bars
   - Check for convergence
   - Compare with experiments

5. Validation:
   - Run multiple independent simulations
   - Check hysteresis (forward/reverse)
   - Compare with other methods
   - Validate with known systems

6. Common pitfalls:
   - Too small s: poor sampling
   - Too large s: loss of accuracy
   - Insufficient sampling
   - Poor initial structures
   - Incorrect topologies
""")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("""
Running EDS from Python:

Quick start:
  from gromos import run_eds
  outputs = run_eds(top, coords, input,
                    num_states=2,
                    smoothness=1.0,
                    steps=100000)

Full control:
  from gromos import EDSSimulation
  eds = EDSSimulation(top, coords, input,
                      num_states=4, smoothness=1.0)
  outputs = eds.run(steps=500000)

Key features:
  ✓ Single simulation, multiple states
  ✓ Better convergence than FEP
  ✓ Free energy differences
  ✓ Python-based analysis (MBAR)
  ✓ Rust performance

Applications:
  - Ligand binding affinities
  - Mutation effects
  - pKa predictions
  - Conformational free energies

Analysis tools:
  - pymbar (MBAR analysis)
  - alchemlyb (unified interface)
  - Custom Python scripts

Next steps:
  - See example 09 for REMD
  - Read EDS paper: J. Chem. Phys. 2002, 116, 8649
  - Install pymbar: pip install pymbar
""")

print("\n" + "=" * 70)
print("Example complete!")
print("=" * 70)
