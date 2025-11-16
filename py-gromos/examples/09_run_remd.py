"""
Example 9: Running REMD (Replica Exchange Molecular Dynamics) from Python
==========================================================================

Replica Exchange MD uses multiple replicas at different temperatures
to enhance sampling of conformational space.
"""

import gromos
import numpy as np

print("=" * 70)
print("Example 9: Running REMD from Python")
print("=" * 70)

# ============================================================================
# What is REMD?
# ============================================================================
print("\n" + "=" * 70)
print("What is REMD?")
print("=" * 70)

print("""
Replica Exchange Molecular Dynamics (REMD) is an enhanced sampling method that:

1. Runs multiple copies (replicas) at different temperatures
2. Periodically attempts to exchange configurations between replicas
3. Overcomes energy barriers at high temperatures
4. Samples low-temperature conformations more efficiently

Exchange criterion (Metropolis):
  P_exchange = min(1, exp(Δβ × ΔE))

Where:
  - Δβ = (1/kT₁ - 1/kT₂)
  - ΔE = (E₁ - E₂)
  - k is Boltzmann constant

Key features:
  - Each replica explores different temperature
  - High-T replicas cross barriers easily
  - Low-T replicas sample target ensemble
  - Exchanges maintain canonical ensemble

Parameters:
  - num_replicas: Number of temperature replicas (8-32 typical)
  - temperatures: Temperature ladder (geometric spacing)
  - exchange_frequency: How often to attempt exchanges (100-1000 steps)
  - steps_per_replica: Simulation length per replica

Applications:
  - Protein folding
  - Conformational sampling
  - Binding free energies
  - Phase transitions
  - Aggregation studies
""")

# ============================================================================
# Method 1: Quick REMD run
# ============================================================================
print("\n" + "=" * 70)
print("Method 1: Quick REMD Simulation")
print("=" * 70)

print("""
from gromos import run_remd

# Run REMD simulation
outputs = run_remd(
    topology="system.top",
    coordinates="equilibrated.cnf",
    input_file="remd.imd",
    temperatures=[300, 310, 320, 330, 340, 350, 360, 370],  # K
    steps=500000,
    exchange_frequency=1000,
    output_prefix="remd"
)

print(f"Replica trajectories: {outputs['trajectories']}")
print(f"Exchange log: {outputs['exchange_log']}")
print(f"Acceptance rates: {outputs['acceptance_rates']}")
print(f"Demultiplexed: {outputs['demultiplexed']}")
""")

# ============================================================================
# Method 2: Using REMD class
# ============================================================================
print("\n" + "=" * 70)
print("Method 2: REMD Class with Custom Parameters")
print("=" * 70)

print("""
from gromos import REMDSimulation

# Create REMD simulation
remd = REMDSimulation(
    topology_file="system.top",
    coordinate_file="start.cnf",
    input_file="remd.imd",
    output_prefix="remd_production",
    temperatures=[300, 310, 320, 330, 340, 350, 360, 370],
    exchange_frequency=1000
)

# Run REMD
outputs = remd.run(
    steps=1000000,
    verbose=True
)

# Analyze exchange statistics
stats = outputs['exchange_statistics']
for i, T in enumerate(remd.temperatures):
    print(f"Replica {i} (T={T}K): acceptance = {stats[i]['acceptance']:.2%}")
""")

# ============================================================================
# Temperature selection
# ============================================================================
print("\n" + "=" * 70)
print("Selecting REMD Temperatures")
print("=" * 70)

print("""
Choosing temperature ladder:
---------------------------

Goal: Achieve 20-40% exchange acceptance between neighbors

Geometric spacing (typical):
  T_i = T_min × (T_max/T_min)^(i/N)

  Where:
  - T_min: Target temperature (e.g., 300 K)
  - T_max: Maximum temperature (e.g., 400-500 K)
  - N: Number of replicas
  - i: Replica index (0 to N-1)

Example for 8 replicas (300-400 K):
  T = [300, 313, 327, 341, 356, 371, 387, 400] K

Optimal spacing depends on:
1. System size (larger = closer spacing needed)
2. Heat capacity (affects ΔE fluctuations)
3. Desired acceptance rate (20-40% is good)

Rule of thumb:
--------------
Temperature increment:
  ΔT ≈ T × sqrt(2 / (3N_dof))

  Where N_dof = degrees of freedom

  Example: 1000 atoms in water
  N_dof ≈ 3 × 1000 = 3000
  ΔT ≈ 300 × sqrt(2/3000) ≈ 7.7 K

  → Use 8-10 K spacing

Number of replicas:
-------------------
Minimum: 4 replicas (barely sufficient)
Typical: 8-16 replicas (good coverage)
Large systems: 32-64 replicas (expensive!)

Cost considerations:
  - Computational cost scales linearly with replicas
  - Need to balance coverage vs. cost
  - 8-16 replicas is sweet spot for most systems

Exchange frequency:
------------------
Too frequent (< 100 steps):
  • Wasted computational effort
  • Overhead of exchange attempts
  • Poor statistics per attempt

Too rare (> 5000 steps):
  • Replicas don't mix well
  • Slow convergence
  • Inefficient sampling

Optimal: 500-1000 steps
  • Good statistics per exchange
  • Reasonable overhead
  • Efficient mixing
""")

# ============================================================================
# Complete REMD workflow
# ============================================================================
print("\n" + "=" * 70)
print("Complete REMD Workflow")
print("=" * 70)

print("""
def run_remd_workflow(topology, coordinates):
    '''
    Complete REMD workflow:
    1. Generate temperature ladder
    2. Equilibrate each replica
    3. Run REMD with exchanges
    4. Demultiplex trajectories
    5. Analyze convergence
    '''
    from gromos import MDSimulation, REMDSimulation
    import numpy as np

    # Step 1: Generate temperature ladder
    print("Step 1: Generating temperature ladder...")
    T_min = 300  # K
    T_max = 400  # K
    num_replicas = 8

    # Geometric spacing
    temperatures = [
        T_min * (T_max / T_min) ** (i / (num_replicas - 1))
        for i in range(num_replicas)
    ]
    print(f"Temperatures: {[f'{T:.1f}' for T in temperatures]} K")

    # Step 2: Equilibrate each replica
    print("Step 2: Equilibrating replicas...")
    equilibrated_configs = []

    for i, T in enumerate(temperatures):
        print(f"  Equilibrating replica {i} at {T:.1f} K...")
        md = MDSimulation(
            topology_file=topology,
            coordinate_file=coordinates,
            input_file=f"equilibrate_T{T:.0f}.imd",  # Temperature-specific input
            output_prefix=f"replica_{i}_eq"
        )
        eq_out = md.run(steps=50000)
        equilibrated_configs.append(eq_out['final_config'])

    # Step 3: Run REMD
    print("Step 3: Running REMD...")
    remd = REMDSimulation(
        topology_file=topology,
        coordinate_file=equilibrated_configs[0],  # Start from equilibrated
        input_file="remd.imd",
        output_prefix="remd",
        temperatures=temperatures,
        exchange_frequency=1000
    )

    remd_out = remd.run(steps=1000000)

    # Step 4: Demultiplex trajectories
    print("Step 4: Demultiplexing trajectories...")
    demux = demultiplex_trajectories(
        remd_out['trajectories'],
        remd_out['exchange_log']
    )

    # Step 5: Analyze convergence
    print("Step 5: Analyzing convergence...")
    stats = analyze_remd(remd_out)

    print(f"\\nREMD Statistics:")
    print(f"  Overall acceptance: {stats['mean_acceptance']:.2%}")
    print(f"  Round-trip time: {stats['round_trip_time']:.1f} ns")
    print(f"  Mixing efficiency: {stats['mixing_efficiency']:.2%}")

    return remd_out, stats

def demultiplex_trajectories(trajectories, exchange_log):
    '''
    Demultiplex: separate by replica index (not temperature)

    Each replica walks through temperature space:
    Time:     0    100   200   300   400   500
    Replica 0: T0 → T1 → T0 → T2 → T0 → T1
    Replica 1: T1 → T0 → T2 → T1 → T2 → T0
    ...

    Demultiplexing extracts continuous trajectory for each replica.
    '''
    import pandas as pd

    # Read exchange log
    exchanges = pd.read_csv(exchange_log, delim_whitespace=True)

    # Track which configuration is in which replica
    # (implementation details depend on GROMOS-RS output format)

    demuxed_trajs = []
    for replica_id in range(len(trajectories)):
        # Extract frames where this replica was sampled
        frames = extract_replica_frames(trajectories, exchanges, replica_id)
        demuxed_trajs.append(frames)

    return demuxed_trajs

def analyze_remd(remd_output):
    '''Analyze REMD performance'''
    import pandas as pd
    import numpy as np

    # Read exchange log
    log = pd.read_csv(remd_output['exchange_log'], delim_whitespace=True)

    # Calculate acceptance rates between neighboring replicas
    acceptances = []
    for i in range(len(remd_output['temperatures']) - 1):
        attempts = log[f'attempts_{i}_{i+1}'].sum()
        accepts = log[f'accepts_{i}_{i+1}'].sum()
        acceptances.append(accepts / attempts if attempts > 0 else 0)

    # Calculate round-trip time (how long for replica to visit all temps)
    round_trip = calculate_round_trip_time(log)

    # Calculate mixing efficiency
    mixing = calculate_mixing_efficiency(log)

    return {
        'mean_acceptance': np.mean(acceptances),
        'acceptances': acceptances,
        'round_trip_time': round_trip,
        'mixing_efficiency': mixing
    }

# Run workflow
# remd_out, stats = run_remd_workflow('system.top', 'start.cnf')
""")

# ============================================================================
# Exchange statistics and analysis
# ============================================================================
print("\n" + "=" * 70)
print("Exchange Statistics and Analysis")
print("=" * 70)

print("""
Key metrics to monitor:

1. Acceptance Rate (per replica pair):
   -----------------------------------
   acceptance = N_accepted / N_attempted

   Optimal: 20-40%
   Too high (>60%): Temperatures too close (inefficient)
   Too low (<10%): Temperatures too far (poor mixing)

   Example output:
     Pair (T=300-310): 35.2% acceptance ✓
     Pair (T=310-320): 32.8% acceptance ✓
     Pair (T=320-330): 28.5% acceptance ✓
     ...

2. Round-Trip Time:
   ----------------
   Time for replica to visit all temperatures and return

   Shorter = better mixing

   Calculation:
     - Track replica 0 over time
     - Measure time to visit T_max and return to T_min
     - Average over multiple round trips

   Typical: 10-100 ns (system dependent)

3. Replica Flow:
   -------------
   Visualize replica movement through temperature space:

   Time (ns)  |  0   1   2   3   4   5   6   7   8
   -----------+-----------------------------------
   Replica 0  |  0 → 1 → 0 → 2 → 1 → 0 → 3 → 2 → 1
   Replica 1  |  1 → 0 → 2 → 0 → 2 → 3 → 1 → 0 → 2
   Replica 2  |  2 → 3 → 1 → 3 → 0 → 1 → 2 → 3 → 0
   ...

   Good flow: Replicas diffuse through temperature space
   Bad flow: Replicas stuck at certain temperatures

4. Potential Energy Distributions:
   -------------------------------
   Plot energy distributions at each temperature:

   Should overlap between neighboring replicas (20-40% overlap)

   Good overlap:        Poor overlap:
     T1  T2               T1      T2
      ╱╲╱╲                ╱╲      ╱╲
     ╱  ╳  ╲              ╱  ╲    ╱  ╲
    ╱  ╱ ╲  ╲            ╱    ╲  ╱    ╲

   Use to diagnose problems and adjust temperature ladder.

5. Mixing Efficiency:
   ------------------
   mixing = σ²_observed / σ²_ideal

   Where:
   - σ²_observed: variance of replica temperature indices
   - σ²_ideal: variance for perfect mixing

   Ideal = 1.0 (perfect mixing)
   Good > 0.7
   Poor < 0.5

Analysis example:
-----------------

def analyze_remd_detailed(remd_output):
    '''Comprehensive REMD analysis'''
    import matplotlib.pyplot as plt
    import numpy as np

    # 1. Plot acceptance rates
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Subplot 1: Acceptance vs replica pair
    axes[0, 0].bar(range(len(acceptances)), acceptances)
    axes[0, 0].axhline(0.2, color='r', linestyle='--', label='Min')
    axes[0, 0].axhline(0.4, color='r', linestyle='--', label='Max')
    axes[0, 0].set_xlabel('Replica Pair')
    axes[0, 0].set_ylabel('Acceptance Rate')
    axes[0, 0].set_title('Exchange Acceptance Rates')
    axes[0, 0].legend()

    # Subplot 2: Replica flow (temperature vs time)
    replica_temps = extract_replica_temperatures(exchange_log)
    for i in range(num_replicas):
        axes[0, 1].plot(replica_temps[:, i], alpha=0.7, label=f'Rep {i}')
    axes[0, 1].set_xlabel('Exchange Attempt')
    axes[0, 1].set_ylabel('Temperature Index')
    axes[0, 1].set_title('Replica Temperature Flow')
    axes[0, 1].legend()

    # Subplot 3: Potential energy distributions
    for i, T in enumerate(temperatures):
        energies = extract_potential_energies(trajectories[i])
        axes[1, 0].hist(energies, bins=50, alpha=0.5, label=f'{T:.0f} K')
    axes[1, 0].set_xlabel('Potential Energy (kJ/mol)')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Energy Distributions')
    axes[1, 0].legend()

    # Subplot 4: Round-trip time distribution
    round_trips = calculate_all_round_trips(exchange_log)
    axes[1, 1].hist(round_trips, bins=20)
    axes[1, 1].set_xlabel('Round-Trip Time (ns)')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Round-Trip Time Distribution')

    plt.tight_layout()
    plt.savefig('remd_analysis.png', dpi=300)

    print("REMD analysis plot saved to remd_analysis.png")

# Run detailed analysis
# analyze_remd_detailed(remd_out)
""")

# ============================================================================
# REMD vs other enhanced sampling methods
# ============================================================================
print("\n" + "=" * 70)
print("REMD vs Other Enhanced Sampling Methods")
print("=" * 70)

print("""
Comparison:

Standard MD:
  Pros:
  + Simple setup
  + Single simulation
  + Easy analysis

  Cons:
  - Slow barrier crossing
  - Poor sampling
  - May get trapped in local minima
  - Long time scales needed

GaMD (Gaussian Accelerated MD):
  Pros:
  + Single simulation
  + Smooth boost potential
  + Can recover original ensemble

  Cons:
  - Parameter tuning required
  - Reweighting needed
  - May over-smooth landscape
  - Numerical instabilities possible

EDS (Enveloping Distribution Sampling):
  Pros:
  + Multiple states in one simulation
  + Good for free energies
  + Better convergence than FEP

  Cons:
  - Limited to predefined states
  - Parameter dependent (smoothness)
  - Requires MBAR analysis

REMD (Replica Exchange MD):
  Pros:
  + No biasing potentials
  + Samples original ensemble at each T
  + Robust and reliable
  + Well-established method
  + Good for folding/unfolding

  Cons:
  - Computationally expensive (N replicas)
  - Need many replicas for large systems
  - Scaling issues for very large systems
  - Exchange overhead

When to use REMD:
  ✓ Protein folding studies
  ✓ Conformational sampling
  ✓ Systems with multiple metastable states
  ✓ When you need unbiased sampling at target T
  ✓ Sufficient computational resources available
  ✗ Very large systems (>100k atoms)
  ✗ Limited computational budget
  ✗ Single specific pathway needed
""")

# ============================================================================
# Advanced REMD variants
# ============================================================================
print("\n" + "=" * 70)
print("Advanced REMD Variants")
print("=" * 70)

print("""
1. Hamiltonian REMD (H-REMD):
   --------------------------
   Exchange between different Hamiltonians instead of temperatures

   Use cases:
   - Free energy calculations
   - Alchemical transformations
   - Partial system scaling

   Example:
     Replica 0: Full interactions
     Replica 1: 90% LJ interactions
     Replica 2: 80% LJ interactions
     ...

2. Solute Tempering (REST):
   -------------------------
   Scale only solute interactions, keep solvent at target T

   Advantages:
   - Fewer replicas needed
   - Solvent remains at biological T
   - More efficient exchange

   Temperature scaling:
     E_total = E_solute/λ + E_solvent + E_solute-solvent/sqrt(λ)

   Where λ = T_0/T_replica

3. Replica Exchange with Solute Scaling (RESS):
   --------------------------------------------
   Combine REST with Hamiltonian scaling

   Scale both temperature and interactions
   Maximum efficiency for protein systems

4. Multi-dimensional REMD:
   -----------------------
   Exchange along multiple parameters simultaneously

   Example: Temperature + pH
     Replica (300K, pH7.0)
     Replica (310K, pH7.0)
     Replica (300K, pH7.5)
     Replica (310K, pH7.5)
     ...

   Cost: Exponential scaling (N_T × N_pH × ...)

5. Adaptive REMD:
   --------------
   Adjust temperature ladder during simulation

   Monitor acceptance rates → adjust temperatures
   Goal: maintain 20-40% acceptance automatically
""")

# ============================================================================
# Best practices
# ============================================================================
print("\n" + "=" * 70)
print("REMD Best Practices")
print("=" * 70)

print("""
1. System preparation:
   - Equilibrate well before REMD
   - Each replica should start from equilibrated structure
   - Check energy conservation
   - Verify temperature control

2. Temperature ladder:
   - Start with geometric spacing
   - Aim for 20-40% acceptance
   - Adjust based on test runs
   - Consider heat capacity of system

3. Simulation length:
   - Minimum: 10 ns per replica (very small systems)
   - Typical: 50-200 ns per replica
   - Goal: Multiple round trips observed
   - Check convergence of properties

4. Exchange frequency:
   - Typical: 500-1000 steps (1-2 ps)
   - Not too frequent (overhead)
   - Not too rare (poor mixing)
   - Balance statistics vs. efficiency

5. Number of replicas:
   - Small peptides: 8-16 replicas
   - Medium proteins: 16-32 replicas
   - Large proteins: 32-64 replicas
   - Cost = N_replicas × single simulation

6. Analysis:
   - Monitor acceptance rates throughout
   - Check replica flow (mixing)
   - Verify round-trip times
   - Analyze energy distributions
   - Use only T_target replica for final analysis

7. Common pitfalls:
   - Too few replicas: poor temperature coverage
   - Too many replicas: wasted resources
   - Wrong spacing: poor acceptance
   - Too short: no round trips observed
   - Not equilibrated: artifacts in early exchanges
   - Forgetting to demultiplex: wrong ensemble

8. Validation:
   - Compare with long standard MD (if possible)
   - Check convergence from different starting structures
   - Verify detailed balance (forward = reverse)
   - Compare with experimental data
""")

# ============================================================================
# Practical example
# ============================================================================
print("\n" + "=" * 70)
print("Practical Example: Protein Folding")
print("=" * 70)

print("""
Case study: Small protein (50 residues) folding

System:
  - 50 residues
  - Explicit solvent
  - ~15,000 atoms total
  - Target: 300 K

REMD setup:
-----------

# 1. Generate temperature ladder
temperatures = generate_geometric_ladder(
    T_min=300,
    T_max=400,
    num_replicas=16
)
# Result: [300, 308, 316, 325, 333, 342, 351, 361, 370, 380, 390, 400] K

# 2. Estimate acceptance
system_size = 15000
dof = 3 * system_size
delta_T = estimate_spacing(T_min, dof)
# delta_T ≈ 8 K (matches our ladder!)

# 3. Setup REMD
from gromos import REMDSimulation

remd = REMDSimulation(
    topology_file="protein_solvated.top",
    coordinate_file="protein_equilibrated.cnf",
    input_file="remd.imd",
    output_prefix="folding_remd",
    temperatures=temperatures,
    exchange_frequency=1000  # Every 1 ps
)

# 4. Run production
outputs = remd.run(
    steps=50_000_000,  # 100 ns per replica
    verbose=True
)

# 5. Analyze
from gromos import analyze_remd

stats = analyze_remd(outputs)
print(f"Mean acceptance: {stats['mean_acceptance']:.2%}")
print(f"Round-trip time: {stats['round_trip_time']:.1f} ns")

# 6. Extract T=300K ensemble
demuxed = demultiplex_trajectories(
    outputs['trajectories'],
    outputs['exchange_log']
)

# Replica 0 trajectory at T=300K
traj_300K = demuxed[0]

# 7. Analyze folding
folding_analysis(traj_300K)

Expected results:
-----------------
- Acceptance: 25-35% (good!)
- Round trips: 5-10 trips in 100 ns
- Folded structures sampled at T=300K
- Unfolded structures sampled at T>350K
- Reversible folding/unfolding observed

Computational cost:
------------------
16 replicas × 100 ns = 1,600 ns total simulation time
On 256 cores: ~2-3 days
""")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("""
Running REMD from Python:

Quick start:
  from gromos import run_remd
  outputs = run_remd(
      topology, coordinates, input_file,
      temperatures=[300, 310, 320, 330, 340, 350],
      steps=1000000,
      exchange_frequency=1000
  )

Full control:
  from gromos import REMDSimulation
  remd = REMDSimulation(
      topology, coordinates, input_file,
      temperatures=temp_ladder,
      exchange_frequency=1000
  )
  outputs = remd.run(steps=1000000)

Key features:
  ✓ Multiple replicas in parallel
  ✓ Automatic exchange attempts
  ✓ Detailed statistics logging
  ✓ Demultiplexing utilities
  ✓ Analysis tools included

Critical parameters:
  - Temperature ladder (geometric spacing)
  - Number of replicas (8-32)
  - Exchange frequency (500-1000 steps)
  - Simulation length (50-200 ns per replica)

Success criteria:
  ✓ Acceptance: 20-40%
  ✓ Multiple round trips observed
  ✓ Good replica mixing
  ✓ Converged properties at T_target

Applications:
  - Protein folding
  - Conformational transitions
  - Ligand binding pathways
  - Aggregation mechanisms
  - Phase transitions

Analysis:
  - Exchange statistics
  - Replica flow diagrams
  - Energy distribution overlap
  - Round-trip time analysis
  - Convergence checks

Next steps:
  - Read REMD papers for theory
  - Test temperature ladder on your system
  - Monitor acceptance rates
  - Adjust parameters as needed
  - Analyze T_target ensemble only
""")

print("\n" + "=" * 70)
print("Example complete!")
print("=" * 70)
