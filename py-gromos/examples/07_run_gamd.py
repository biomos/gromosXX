"""
Example 7: Running GaMD (Gaussian Accelerated MD) from Python
==============================================================

Gaussian Accelerated Molecular Dynamics adds a harmonic boost potential
to smooth the energy landscape, accelerating sampling of rare events.
"""

import gromos
import numpy as np

print("=" * 70)
print("Example 7: Running GaMD from Python")
print("=" * 70)

# ============================================================================
# What is GaMD?
# ============================================================================
print("\n" + "=" * 70)
print("What is GaMD?")
print("=" * 70)

print("""
Gaussian Accelerated Molecular Dynamics (GaMD) is an enhanced sampling
method that:

1. Adds a harmonic boost potential to the system
2. Smooths the energy landscape
3. Accelerates rare events (transitions, folding, binding)
4. Enables reweighting to canonical ensemble

Boost potential:
  ΔV(r) = 1/2 * k * (E - V(r))²  when V(r) < E
  ΔV(r) = 0                       when V(r) ≥ E

Where:
  - E is the threshold energy
  - k is the boost force constant
  - V(r) is the potential energy

Key parameter:
  - sigma0: controls boost magnitude (6-8 kJ/mol typical)

Applications:
  - Protein folding/unfolding
  - Ligand binding/unbinding
  - Conformational transitions
  - Allosteric mechanisms
""")

# ============================================================================
# Method 1: Quick GaMD run
# ============================================================================
print("\n" + "=" * 70)
print("Method 1: Quick GaMD Simulation")
print("=" * 70)

print("""
from gromos import run_gamd

# Run GaMD with default parameters
outputs = run_gamd(
    topology="protein.top",
    coordinates="equilibrated.cnf",
    input_file="gamd.imd",
    equilibration_steps=10000,   # Collect energy statistics
    production_steps=100000,      # Run with boost
    sigma0=6.0,                   # Boost parameter (kJ/mol)
    threshold_mode='lower',       # Apply boost when E < threshold
    output_prefix="gamd"
)

print(f"Trajectory: {outputs['trajectory']}")
print(f"Energy: {outputs['energy']}")
print(f"Boost potential: {outputs['boost']}")
print(f"Statistics: {outputs['statistics']}")
""")

# ============================================================================
# Method 2: Using GaMD class
# ============================================================================
print("\n" + "=" * 70)
print("Method 2: GaMD Class with Custom Parameters")
print("=" * 70)

print("""
from gromos import GaMDSimulation

# Create GaMD simulation
gamd = GaMDSimulation(
    topology_file="system.top",
    coordinate_file="start.cnf",
    input_file="gamd.imd",
    output_prefix="gamd_production",
    sigma0=6.0,           # Standard deviation of boost (kJ/mol)
    threshold_mode='lower'  # or 'upper'
)

# Run GaMD
outputs = gamd.run(
    equilibration_steps=50000,  # Longer equilibration
    production_steps=500000,     # Longer production
    verbose=True
)

# Outputs include boost potential trajectory
boost_file = outputs['boost']
stats_file = outputs['statistics']
""")

# ============================================================================
# Parameter selection
# ============================================================================
print("\n" + "=" * 70)
print("Selecting GaMD Parameters")
print("=" * 70)

print("""
Choosing sigma0:
----------------
- Smaller sigma0 (4-5 kJ/mol): Gentle acceleration, more accurate
- Larger sigma0 (7-8 kJ/mol): Stronger acceleration, faster sampling
- Typical: 6 kJ/mol (good balance)

Rule of thumb:
  sigma0 ≈ 0.2 * sqrt(N_atoms) kJ/mol

For different system sizes:
  - 100 atoms:    sigma0 ≈ 2 kJ/mol
  - 1,000 atoms:  sigma0 ≈ 6 kJ/mol
  - 10,000 atoms: sigma0 ≈ 20 kJ/mol

Threshold mode:
---------------
- 'lower': Apply boost when V(r) < E (typical)
  → Fills energy wells, accelerates escapes

- 'upper': Apply boost when V(r) > E (less common)
  → Reduces barriers

Equilibration:
--------------
Minimum steps = 10 * autocorrelation_time
Typical: 10,000 - 100,000 steps

Production:
-----------
Depends on system and timescale of interest
Typical: 100,000 - 10,000,000 steps
""")

# ============================================================================
# Complete GaMD workflow
# ============================================================================
print("\n" + "=" * 70)
print("Complete GaMD Workflow")
print("=" * 70)

print("""
def run_gamd_workflow(topology, coordinates):
    '''
    Complete GaMD workflow:
    1. Conventional MD equilibration
    2. GaMD equilibration (collect statistics)
    3. GaMD production
    4. Reweighting analysis
    '''
    from gromos import MDSimulation, GaMDSimulation
    import numpy as np

    # Step 1: Conventional MD equilibration
    print("Step 1: Conventional MD equilibration...")
    md_eq = MDSimulation(
        topology_file=topology,
        coordinate_file=coordinates,
        input_file="md_equilibration.imd",
        output_prefix="md_eq"
    )
    md_out = md_eq.run(steps=50000)

    # Step 2-3: GaMD (equilibration + production combined)
    print("Step 2-3: GaMD simulation...")
    gamd = GaMDSimulation(
        topology_file=topology,
        coordinate_file=str(md_out['final_config']),
        input_file="gamd.imd",
        output_prefix="gamd",
        sigma0=6.0,
        threshold_mode='lower'
    )
    gamd_out = gamd.run(
        equilibration_steps=50000,
        production_steps=500000
    )

    # Step 4: Analyze and reweight
    print("Step 4: Analyzing GaMD results...")
    analyze_gamd_output(gamd_out)

    return gamd_out

def analyze_gamd_output(outputs):
    '''Analyze GaMD simulation'''
    import pandas as pd
    import matplotlib.pyplot as plt

    # Read boost potential
    boost_data = pd.read_csv(
        outputs['boost'],
        delim_whitespace=True,
        names=['time', 'boost', 'potential']
    )

    # Plot boost potential
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    axes[0].plot(boost_data['time'], boost_data['boost'])
    axes[0].set_xlabel('Time (ps)')
    axes[0].set_ylabel('Boost Potential (kJ/mol)')
    axes[0].set_title('GaMD Boost Potential')

    axes[1].hist(boost_data['boost'], bins=50, alpha=0.7)
    axes[1].set_xlabel('Boost Potential (kJ/mol)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Boost Distribution')

    plt.tight_layout()
    plt.savefig('gamd_boost_analysis.png', dpi=300)

    print(f"Average boost: {boost_data['boost'].mean():.2f} kJ/mol")
    print(f"Std boost: {boost_data['boost'].std():.2f} kJ/mol")

    # Reweighting (simplified)
    # Full reweighting requires cumulant expansion
    weights = np.exp(boost_data['boost'] / (kb * T))
    weights_normalized = weights / weights.sum()

    return weights_normalized

# Run workflow
# gamd_results = run_gamd_workflow('system.top', 'equilibrated.cnf')
""")

# ============================================================================
# Reweighting GaMD results
# ============================================================================
print("\n" + "=" * 70)
print("Reweighting GaMD Trajectories")
print("=" * 70)

print("""
GaMD modifies the potential energy, so trajectories must be reweighted
to recover canonical ensemble properties.

Reweighting formula (cumulant expansion to 2nd order):
  w(t) = exp(β * ΔV(t)) * exp(-β * <ΔV> + β²/2 * σ²ΔV)

Where:
  - β = 1/(kB*T)
  - ΔV(t) = boost potential at time t
  - <ΔV> = average boost
  - σ²ΔV = variance of boost

Implementation:
def reweight_gamd(boost_file, temperature=300):
    '''Reweight GaMD trajectory'''
    import numpy as np
    import pandas as pd

    # Constants
    kb = 0.00831446261815324  # kJ/(mol*K)
    beta = 1.0 / (kb * temperature)

    # Read boost potential
    data = pd.read_csv(boost_file, delim_whitespace=True)
    boost = data['boost'].values

    # Calculate statistics
    avg_boost = boost.mean()
    var_boost = boost.var()

    # Calculate weights (cumulant expansion)
    weights = np.exp(
        beta * boost - beta * avg_boost + 0.5 * beta**2 * var_boost
    )

    # Normalize
    weights /= weights.sum()

    return weights

# Calculate reweighted average
weights = reweight_gamd('gamd_boost.dat', temperature=300)

# Reweighted property
property_values = [...]  # Your property of interest
reweighted_avg = np.average(property_values, weights=weights)

# Reweighted free energy
# ΔG = -kT * ln(P_reweighted)
""")

# ============================================================================
# GaMD vs Standard MD
# ============================================================================
print("\n" + "=" * 70)
print("GaMD vs Standard MD")
print("=" * 70)

print("""
Comparison:

Standard MD:
  Pros:
  + Exact canonical ensemble
  + No reweighting needed
  + Well-established

  Cons:
  - Slow sampling of rare events
  - May not escape energy wells
  - Limited timescales

GaMD:
  Pros:
  + Faster sampling (10-100× speedup)
  + Escapes energy wells
  + No predefined reaction coordinate
  + Enables reweighting

  Cons:
  - Requires reweighting
  - Parameter tuning (sigma0)
  - More complex analysis

When to use GaMD:
  ✓ Rare events (folding, binding)
  ✓ Conformational transitions
  ✓ High barriers
  ✓ Limited simulation time
  ✗ Already fast processes
  ✗ Need exact dynamics
""")

# ============================================================================
# Best practices
# ============================================================================
print("\n" + "=" * 70)
print("GaMD Best Practices")
print("=" * 70)

print("""
1. System preparation:
   - Equilibrate well with standard MD first
   - Ensure system is stable
   - Remove bad contacts

2. Parameter selection:
   - Start with sigma0 = 6.0 kJ/mol
   - Adjust based on system size
   - Use 'lower' threshold (typical)

3. Equilibration:
   - Minimum 10,000 steps
   - Check energy convergence
   - Verify boost distribution

4. Production:
   - Run long enough (100,000+ steps)
   - Monitor boost potential
   - Check for convergence

5. Analysis:
   - Always reweight!
   - Use cumulant expansion
   - Calculate error bars
   - Compare with standard MD

6. Validation:
   - Run multiple independent simulations
   - Check convergence
   - Validate with experimental data
""")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("""
Running GaMD from Python:

Quick start:
  from gromos import run_gamd
  outputs = run_gamd(top, coords, input,
                     equilibration_steps=10000,
                     production_steps=100000,
                     sigma0=6.0)

Full control:
  from gromos import GaMDSimulation
  gamd = GaMDSimulation(top, coords, input, sigma0=6.0)
  outputs = gamd.run(equilibration_steps=50000,
                     production_steps=500000)

Key features:
  ✓ Accelerated sampling (10-100× speedup)
  ✓ No reaction coordinate needed
  ✓ Reweighting to canonical ensemble
  ✓ Python-based analysis
  ✓ Rust performance

Applications:
  - Protein folding
  - Ligand binding
  - Conformational changes
  - Allosteric mechanisms

Next steps:
  - See example 08 for EDS
  - See example 09 for REMD
  - Read the GaMD paper: J. Chem. Theory Comput. 2015, 11, 3584
""")

print("\n" + "=" * 70)
print("Example complete!")
print("=" * 70)
