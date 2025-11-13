"""
Example 4: Advanced Sampling Methods
=====================================

Demonstrates enhanced sampling techniques: GaMD, EDS, and REMD.
"""

import gromos

print("=" * 60)
print("GROMOS-RS: Advanced Sampling Methods")
print("=" * 60)

# ============================================================================
# Gaussian Accelerated Molecular Dynamics (GaMD)
# ============================================================================
print("\n" + "=" * 60)
print("1. Gaussian Accelerated MD (GaMD)")
print("=" * 60)

print("\nWhat is GaMD?")
print("  - Adds harmonic boost potential to smooth energy landscape")
print("  - Accelerates sampling of rare events")
print("  - Enables reweighting to recover canonical ensemble")
print("  - No pre-defined reaction coordinates needed")

print("\n--- GaMD Setup ---")
gamd_params = gromos.GamdParameters(
    sigma0=6.0,  # Standard deviation of energy boost (kJ/mol)
    threshold_mode='lower'  # Apply boost when E < E_threshold
)
print(f"{gamd_params}")

gamd = gromos.GamdRunner(gamd_params)
print(f"{gamd}")

print("\nGaMD Parameters:")
print(f"  sigma0 = {6.0} kJ/mol (controls boost magnitude)")
print(f"  threshold_mode = 'lower' (boost below threshold)")
print("  Alternative: 'upper' (boost above threshold)")

print("\nTypical GaMD Workflow:")
print("  1. Equilibration phase: collect energy statistics")
print("  2. GaMD phase: apply boost potential")
print("  3. Analysis: reweight frames to canonical ensemble")

print("\nApplications:")
print("  - Protein folding/unfolding")
print("  - Ligand binding/unbinding")
print("  - Conformational transitions")
print("  - Allosteric mechanisms")

# ============================================================================
# Enveloping Distribution Sampling (EDS)
# ============================================================================
print("\n" + "=" * 60)
print("2. Enveloping Distribution Sampling (EDS)")
print("=" * 60)

print("\nWhat is EDS?")
print("  - Samples multiple end-states simultaneously")
print("  - Uses smoothed envelope potential")
print("  - Efficient for free energy calculations")
print("  - Better phase space overlap than standard FEP")

print("\n--- EDS Setup ---")
num_states = 4
smoothness = 1.0

eds_params = gromos.EDSParameters(
    num_states=num_states,
    smoothness=smoothness
)
print(f"{eds_params}")

eds = gromos.EDSRunner(eds_params)
print(f"{eds}")

print("\nEDS Parameters:")
print(f"  num_states = {num_states} (number of end-states)")
print(f"  smoothness s = {smoothness} kJ/mol (controls envelope width)")
print("  Larger s → smoother envelope, better sampling")
print("  Smaller s → closer to individual states")

print("\nEDS Envelope Potential:")
print("  V_EDS = -s * ln(Σᵢ exp(-Vᵢ/s))")
print("  where Vᵢ are the end-state potentials")

print("\nApplications:")
print("  - Alchemical free energy calculations")
print("  - Ligand binding affinities")
print("  - Mutation free energies")
print("  - pKa predictions")

print("\nAdvantages over standard FEP:")
print("  + Better sampling of intermediate states")
print("  + Faster convergence")
print("  + Can handle larger perturbations")
print("  + Multiple states in single simulation")

# ============================================================================
# Replica Exchange Molecular Dynamics (REMD)
# ============================================================================
print("\n" + "=" * 60)
print("3. Replica Exchange MD (REMD)")
print("=" * 60)

print("\nWhat is REMD?")
print("  - Runs multiple replicas at different conditions")
print("  - Periodically attempts to exchange configurations")
print("  - Enhances sampling via temperature/parameter random walk")
print("  - Each replica samples canonical ensemble")

print("\n--- REMD Setup ---")
num_replicas = 8
exchange_interval = 1000  # steps

remd = gromos.ReplicaController(
    num_replicas=num_replicas,
    exchange_interval=exchange_interval
)
print(f"{remd}")
print(f"Number of replicas: {remd.num_replicas()}")

print("\nREMD Parameters:")
print(f"  num_replicas = {num_replicas}")
print(f"  exchange_interval = {exchange_interval} steps")
print("  Exchange every {:.1f} ps (at dt=0.002 ps)".format(
    exchange_interval * 0.002
))

print("\nTemperature-REMD Example:")
print("  Replica  Temperature")
print("  " + "-" * 25)
temperatures = [300, 310, 320, 330, 345, 360, 375, 390]
for i, T in enumerate(temperatures):
    print(f"    {i}        {T} K")

print("\nTemperature Spacing:")
print("  - Geometric spacing often used: T(i+1)/T(i) ≈ constant")
print("  - Acceptance ratio target: 20-30%")
print("  - More replicas → better sampling but more expensive")

print("\nExchange Schemes:")
print("  - Nearest neighbor: only adjacent T can exchange")
print("  - Random pairs: any pair can exchange")
print("  - Metropolis criterion: accept if exp(-ΔE/kT) > random")

print("\nApplications:")
print("  - Protein folding")
print("  - Peptide structure prediction")
print("  - Binding free energies")
print("  - Overcoming energy barriers")

# ============================================================================
# Method Comparison
# ============================================================================
print("\n" + "=" * 60)
print("Comparison of Methods")
print("=" * 60)

print("\n{:<12} {:<30} {:<30}".format("Method", "Best For", "Computational Cost"))
print("-" * 72)
methods = [
    ("GaMD", "Rare events, no predefined CV", "1× (single simulation)"),
    ("EDS", "Free energy calculations", "1× (all states together)"),
    ("REMD", "Temperature-dependent sampling", f"{num_replicas}× (multiple replicas)"),
]

for method, best_for, cost in methods:
    print("{:<12} {:<30} {:<30}".format(method, best_for, cost))

print("\n{:<12} {:<20} {:<20}".format("Method", "Pre-requisites", "Post-processing"))
print("-" * 52)
details = [
    ("GaMD", "Energy statistics", "Reweighting"),
    ("EDS", "End-state topologies", "WHAM/MBAR"),
    ("REMD", "Temperature series", "Demultiplexing"),
]

for method, prereq, postproc in details:
    print("{:<12} {:<20} {:<20}".format(method, prereq, postproc))

# ============================================================================
# Choosing the Right Method
# ============================================================================
print("\n" + "=" * 60)
print("Choosing the Right Method")
print("=" * 60)

print("\nUse GaMD when:")
print("  ✓ You need to explore conformational space")
print("  ✓ Transition barriers are high")
print("  ✓ You don't know the reaction coordinate")
print("  ✓ You want to stay at target temperature")
print("  ✗ Don't use for: precise free energies (use EDS/REMD)")

print("\nUse EDS when:")
print("  ✓ Calculating free energy differences")
print("  ✓ Alchemical transformations")
print("  ✓ Multiple end-states to compare")
print("  ✓ You want better convergence than FEP")
print("  ✗ Don't use for: simple conformational sampling (use GaMD)")

print("\nUse REMD when:")
print("  ✓ System has many local minima")
print("  ✓ Temperature-dependent phenomena")
print("  ✓ You have computational resources")
print("  ✓ You need ensemble at each temperature")
print("  ✗ Don't use for: very large systems (too expensive)")

print("\n" + "=" * 60)
print("Advanced sampling examples completed successfully!")
print("=" * 60)
