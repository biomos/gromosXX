"""
Example 10: Running TI (Thermodynamic Integration) from Python
==============================================================

Thermodynamic Integration calculates free energy differences by
integrating the derivative of the Hamiltonian with respect to λ.
"""

import gromos
import numpy as np

print("=" * 70)
print("Example 10: Running TI from Python")
print("=" * 70)

# ============================================================================
# What is TI?
# ============================================================================
print("\n" + "=" * 70)
print("What is TI (Thermodynamic Integration)?")
print("=" * 70)

print("""
Thermodynamic Integration (TI) is a free energy method that:

1. Defines a coupling parameter λ (0 → 1)
2. Runs simulations at multiple λ values
3. Calculates ∂H/∂λ at each λ
4. Integrates to obtain ΔG

Fundamental equation:
  ΔG = ∫₀¹ ⟨∂H/∂λ⟩_λ dλ

Where:
  - H(λ) is the Hamiltonian at coupling parameter λ
  - ⟨∂H/∂λ⟩_λ is the ensemble average at fixed λ
  - Integration is from state A (λ=0) to state B (λ=1)

Hamiltonian interpolation:
  H(λ) = (1-λ) × H_A + λ × H_B

Or more generally (soft-core):
  H(λ) = H_bonded + H_LJ(λ) + H_Coulomb(λ)

Key parameters:
  - num_windows: Number of λ values (11-21 typical)
  - λ_values: Distribution of λ (linear, nonlinear, or optimized)
  - steps_per_window: Simulation length at each λ (50k-500k steps)
  - integration_method: Trapezoid, Simpson's, or Gaussian quadrature

Applications:
  - Absolute binding free energies
  - Relative binding free energies
  - Solvation free energies
  - Mutation effects
  - Conformational free energies
""")

# ============================================================================
# Method 1: Quick TI run
# ============================================================================
print("\n" + "=" * 70)
print("Method 1: Quick TI Simulation")
print("=" * 70)

print("""
from gromos import run_ti

# Run TI simulation with automatic λ selection
outputs = run_ti(
    topology="system.top",
    coordinates="equilibrated.cnf",
    input_file="ti.imd",
    num_windows=11,              # 11 λ values
    steps_per_window=100000,     # 100k steps each
    lambda_spacing="linear",     # Linear spacing 0.0, 0.1, ..., 1.0
    output_prefix="ti"
)

print(f"Lambda values: {outputs['lambda_values']}")
print(f"dH/dλ values: {outputs['dhdl_values']}")
print(f"Free energy: {outputs['free_energy']:.2f} ± {outputs['error']:.2f} kJ/mol")
print(f"Individual window outputs: {outputs['window_outputs']}")
""")

# ============================================================================
# Method 2: Using TI class
# ============================================================================
print("\n" + "=" * 70)
print("Method 2: TI Class with Custom Parameters")
print("=" * 70)

print("""
from gromos import TISimulation
import numpy as np

# Create TI simulation with custom λ values
# Use more points near endpoints (better for soft-core)
lambda_values = np.array([
    0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 0.95, 1.0
])

ti = TISimulation(
    topology_file="system.top",
    coordinate_file="start.cnf",
    input_file="ti.imd",
    output_prefix="ti_production",
    lambda_values=lambda_values
)

# Run all λ-windows
outputs = ti.run(
    steps_per_window=500000,
    equilibration_steps=50000,  # Equilibrate each window first
    verbose=True
)

# Calculate free energy by integration
dG, dG_error = ti.integrate(method="trapezoid")
print(f"ΔG = {dG:.2f} ± {dG_error:.2f} kJ/mol")
""")

# ============================================================================
# Lambda spacing strategies
# ============================================================================
print("\n" + "=" * 70)
print("Selecting λ Values and Spacing")
print("=" * 70)

print("""
Choosing λ spacing:
-------------------

1. Linear spacing (simple, often suboptimal):
   λ = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

   Pros: Easy to understand
   Cons: May undersample regions with large ∂H/∂λ

2. Nonlinear spacing (more points near endpoints):
   For soft-core potentials, concentrate points at λ≈0 and λ≈1

   λ = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]

   Or use polynomial spacing:
   λ_i = (i/N)^p  for creation (p≈2)
   λ_i = 1 - ((N-i)/N)^p  for annihilation

3. Optimized spacing:
   Based on variance of ∂H/∂λ

   More points where σ(∂H/∂λ) is large
   Fewer points where σ(∂H/∂λ) is small

   Requires pilot calculations to determine

4. Adaptive spacing:
   Start with coarse grid
   → Run simulations
   → Identify regions with large fluctuations
   → Add more λ points
   → Rerun

Number of windows:
------------------
Minimum: 9 windows (poor accuracy)
Typical: 11-21 windows (good balance)
Large perturbations: 21-51 windows (expensive but accurate)

Rule of thumb:
  - Simple transformations: 11 windows
  - Charge changes: 15-21 windows
  - Soft-core (LJ): 21-31 windows
  - Electrostatics + LJ: 31-51 windows

Integration methods:
--------------------

1. Trapezoid rule (simplest):
   ΔG = Σᵢ [(∂H/∂λ)ᵢ + (∂H/∂λ)ᵢ₊₁] / 2 × Δλᵢ

2. Simpson's rule (better):
   Requires odd number of points
   ΔG = Δλ/3 × [f₀ + 4f₁ + 2f₂ + 4f₃ + ... + f_n]

3. Gaussian quadrature (most accurate):
   Uses optimal weights
   ΔG = Σᵢ wᵢ × (∂H/∂λ)ᵢ

4. Spline integration:
   Fit cubic spline to ∂H/∂λ vs λ
   Integrate spline analytically

Error estimation:
-----------------
1. Block averaging per window
2. Bootstrapping across windows
3. Propagate uncertainties through integration
4. Compare different integration methods
""")

# ============================================================================
# Complete TI workflow
# ============================================================================
print("\n" + "=" * 70)
print("Complete TI Workflow")
print("=" * 70)

print("""
def run_ti_workflow(topology_A, topology_B, coordinates):
    '''
    Complete TI workflow for free energy calculation:
    1. Prepare λ-dependent topologies
    2. Equilibrate each λ window
    3. Run production at each λ
    4. Calculate ∂H/∂λ for each window
    5. Integrate to get ΔG
    6. Error analysis
    '''
    from gromos import MDSimulation, TISimulation
    import numpy as np
    from scipy import integrate

    # Step 1: Define λ schedule
    print("Step 1: Defining λ schedule...")

    # Use more points near endpoints for soft-core
    lambda_values = np.concatenate([
        np.array([0.0]),
        np.linspace(0.05, 0.45, 9),
        np.linspace(0.55, 0.95, 9),
        np.array([1.0])
    ])
    print(f"Using {len(lambda_values)} λ windows")

    # Step 2: Equilibrate each window
    print("Step 2: Equilibrating λ windows...")

    equilibrated_configs = []
    for i, lam in enumerate(lambda_values):
        print(f"  Window {i}: λ = {lam:.3f}")

        # Create λ-specific input file
        create_lambda_input(lam, f"ti_lambda{lam:.3f}.imd")

        # Equilibrate
        md = MDSimulation(
            topology_file=topology_A,  # Hybrid topology
            coordinate_file=coordinates,
            input_file=f"ti_lambda{lam:.3f}.imd",
            output_prefix=f"ti_eq_lambda{lam:.3f}"
        )
        eq_out = md.run(steps=50000)
        equilibrated_configs.append(eq_out['final_config'])

    # Step 3: Production runs
    print("Step 3: Running production simulations...")

    dhdl_values = []
    dhdl_errors = []

    for i, lam in enumerate(lambda_values):
        print(f"  Window {i}: λ = {lam:.3f}")

        md = MDSimulation(
            topology_file=topology_A,
            coordinate_file=equilibrated_configs[i],
            input_file=f"ti_lambda{lam:.3f}.imd",
            output_prefix=f"ti_prod_lambda{lam:.3f}"
        )
        prod_out = md.run(steps=500000)

        # Extract ∂H/∂λ from output
        dhdl, dhdl_err = extract_dhdl(prod_out['energy'])
        dhdl_values.append(dhdl)
        dhdl_errors.append(dhdl_err)

        print(f"    ⟨∂H/∂λ⟩ = {dhdl:.2f} ± {dhdl_err:.2f} kJ/mol")

    # Step 4: Integrate
    print("Step 4: Integrating ∂H/∂λ...")

    dhdl_values = np.array(dhdl_values)
    dhdl_errors = np.array(dhdl_errors)

    # Method 1: Trapezoid rule
    dG_trap = np.trapz(dhdl_values, lambda_values)

    # Method 2: Simpson's rule (if odd number of points)
    if len(lambda_values) % 2 == 1:
        dG_simp = integrate.simps(dhdl_values, lambda_values)
    else:
        dG_simp = dG_trap

    # Method 3: Cubic spline
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(lambda_values, dhdl_values)
    dG_spline = cs.integrate(0, 1)

    # Error estimation (propagate uncertainties)
    # Simple: assume uncorrelated and use trapezoid
    dG_error = np.sqrt(np.sum((dhdl_errors * np.diff(lambda_values,
                      prepend=0, append=0) / 2)**2))

    # Step 5: Report results
    print(f"\\nFree Energy Results:")
    print(f"  Trapezoid:  ΔG = {dG_trap:.2f} ± {dG_error:.2f} kJ/mol")
    print(f"  Simpson:    ΔG = {dG_simp:.2f} kJ/mol")
    print(f"  Spline:     ΔG = {dG_spline:.2f} kJ/mol")
    print(f"  Mean:       ΔG = {np.mean([dG_trap, dG_simp, dG_spline]):.2f} kJ/mol")

    return {
        'lambda_values': lambda_values,
        'dhdl_values': dhdl_values,
        'dhdl_errors': dhdl_errors,
        'free_energy': dG_trap,
        'free_energy_error': dG_error
    }

def extract_dhdl(energy_file):
    '''Extract ∂H/∂λ from energy file'''
    import pandas as pd

    # Read energy file
    data = pd.read_csv(energy_file, delim_whitespace=True, comment='#')

    # Extract ∂H/∂λ column
    dhdl_series = data['dHdl']

    # Discard equilibration (first 20%)
    n_equil = len(dhdl_series) // 5
    dhdl_series = dhdl_series[n_equil:]

    # Calculate mean and error (block averaging)
    dhdl_mean = dhdl_series.mean()
    dhdl_error = block_average_error(dhdl_series)

    return dhdl_mean, dhdl_error

def block_average_error(data, n_blocks=5):
    '''Calculate error using block averaging'''
    import numpy as np

    block_size = len(data) // n_blocks
    block_means = []

    for i in range(n_blocks):
        block = data[i*block_size:(i+1)*block_size]
        block_means.append(block.mean())

    # Standard error of block means
    return np.std(block_means) / np.sqrt(n_blocks)

def create_lambda_input(lam, output_file):
    '''Create input file for specific λ value'''
    # Modify PERTURBATION block in input file
    # Set RLAM = lam
    # (implementation depends on GROMOS input format)
    pass

# Run workflow
# results = run_ti_workflow('hybrid.top', 'stateB.top', 'start.cnf')
""")

# ============================================================================
# TI curve analysis
# ============================================================================
print("\n" + "=" * 70)
print("TI Curve Analysis")
print("=" * 70)

print("""
Analyzing ∂H/∂λ vs λ curve:
---------------------------

The shape reveals important information:

1. Smooth curve (ideal):
   ∂H/∂λ
    │     ╱‾‾‾‾╲
    │    ╱       ╲
    │   ╱         ╲
    │  ╱           ╲
    └──────────────────→ λ
   0.0             1.0

   → Good sampling
   → Sufficient windows
   → Accurate integration

2. Noisy curve (poor sampling):
   ∂H/∂λ
    │    ╱╲    ╱╲╱╲
    │   ╱  ╲  ╱    ╲╱╲
    │  ╱    ╲╱
    │ ╱
    └──────────────────→ λ

   → Increase simulation time per window
   → Check for phase transitions
   → Verify equilibration

3. Sharp peaks (need more windows):
   ∂H/∂λ
    │      │
    │      │
    │   ╱‾‾│‾‾╲
    │  ╱   │   ╲
    └──────┼──────→ λ
         0.5

   → Add more λ points near peak
   → Use soft-core potentials
   → Check for singularities

4. Hysteresis (non-equilibrium):
   Forward and reverse don't match
   → Increase equilibration
   → Increase production time
   → Check for slow degrees of freedom

Visualization and diagnostics:
------------------------------

import matplotlib.pyplot as plt
import numpy as np

def plot_ti_analysis(lambda_values, dhdl_values, dhdl_errors):
    '''Comprehensive TI curve visualization'''

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Subplot 1: ∂H/∂λ vs λ with error bars
    axes[0, 0].errorbar(lambda_values, dhdl_values, yerr=dhdl_errors,
                        fmt='o-', capsize=5, label='Data')

    # Fit cubic spline for smooth curve
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(lambda_values, dhdl_values)
    lam_fine = np.linspace(0, 1, 200)
    axes[0, 0].plot(lam_fine, cs(lam_fine), '--', alpha=0.5, label='Spline fit')

    axes[0, 0].set_xlabel('λ')
    axes[0, 0].set_ylabel('⟨∂H/∂λ⟩ (kJ/mol)')
    axes[0, 0].set_title('TI Curve')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # Subplot 2: Cumulative integral (convergence)
    dG_cumulative = np.zeros(len(lambda_values))
    for i in range(1, len(lambda_values)):
        dG_cumulative[i] = dG_cumulative[i-1] + \
            0.5 * (dhdl_values[i] + dhdl_values[i-1]) * \
            (lambda_values[i] - lambda_values[i-1])

    axes[0, 1].plot(lambda_values, dG_cumulative, 'o-')
    axes[0, 1].axhline(dG_cumulative[-1], color='r', linestyle='--',
                       label=f'Final ΔG = {dG_cumulative[-1]:.2f} kJ/mol')
    axes[0, 1].set_xlabel('λ')
    axes[0, 1].set_ylabel('Cumulative ΔG (kJ/mol)')
    axes[0, 1].set_title('Integration Convergence')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # Subplot 3: Error bars magnitude
    axes[1, 0].bar(lambda_values, dhdl_errors)
    axes[1, 0].set_xlabel('λ')
    axes[1, 0].set_ylabel('Error (kJ/mol)')
    axes[1, 0].set_title('Statistical Uncertainty per Window')
    axes[1, 0].grid(True, alpha=0.3)

    # Subplot 4: Variance of ∂H/∂λ (from errors)
    # High variance → need more sampling or more windows
    relative_error = dhdl_errors / np.abs(dhdl_values)
    axes[1, 1].bar(lambda_values, relative_error * 100)
    axes[1, 1].axhline(10, color='r', linestyle='--', label='10% threshold')
    axes[1, 1].set_xlabel('λ')
    axes[1, 1].set_ylabel('Relative Error (%)')
    axes[1, 1].set_title('Relative Uncertainty per Window')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('ti_analysis.png', dpi=300)
    print("TI analysis plot saved to ti_analysis.png")

# Usage:
# plot_ti_analysis(results['lambda_values'],
#                  results['dhdl_values'],
#                  results['dhdl_errors'])

Convergence checks:
-------------------

1. Time convergence (per window):
   Plot ⟨∂H/∂λ⟩ vs simulation time
   Should plateau after equilibration

2. Forward vs reverse:
   Run TI in both directions (0→1 and 1→0)
   Should give same result (within error)

3. Number of windows:
   Repeat with more windows
   ΔG should converge

4. Integration method:
   Compare trapezoid, Simpson, spline
   Should agree within error bars
""")

# ============================================================================
# TI vs other free energy methods
# ============================================================================
print("\n" + "=" * 70)
print("TI vs Other Free Energy Methods")
print("=" * 70)

print("""
Comparison:

FEP (Free Energy Perturbation):
  Equation: ΔG = -kT ln⟨exp(-ΔU/kT)⟩

  Pros:
  + Can use fewer windows
  + Direct exponential averaging

  Cons:
  - Poor phase space overlap → large errors
  - Sensitive to outliers
  - Difficult for large perturbations
  - Slow convergence

TI (Thermodynamic Integration):
  Equation: ΔG = ∫₀¹ ⟨∂H/∂λ⟩ dλ

  Pros:
  + Smooth free energy profile
  + More stable than FEP
  + Good error estimates
  + Less sensitive to outliers
  + Well-established theory

  Cons:
  - Requires many windows (11-51)
  - Numerical integration errors
  - Need good λ spacing
  - Expensive (many simulations)

BAR (Bennett Acceptance Ratio):
  Equation: Uses forward and reverse work

  Pros:
  + Most efficient use of data
  + Fewer windows than TI
  + Good error estimates

  Cons:
  - More complex analysis
  - Need overlap between windows
  - Requires pymbar or similar

EDS (Enveloping Distribution Sampling):
  Equation: V_EDS = -s ln(Σ exp(-V_i/s))

  Pros:
  + Single simulation
  + Multiple end-states
  + Better overlap

  Cons:
  - Parameter tuning (smoothness)
  - Requires MBAR analysis
  - Limited to discrete states

When to use TI:
  ✓ Smooth transformations needed
  ✓ Want detailed free energy profile
  ✓ Need robust error estimates
  ✓ Standard method for comparison
  ✓ Well-documented protocol
  ✗ Limited computational resources
  ✗ Need only endpoint free energy difference
  ✗ Simple two-state system (use BAR instead)

Hybrid approaches:
------------------
Often combine methods:
- TI for difficult regions (endpoints)
- BAR for smooth regions (middle)
- MBAR for optimal data usage
- Use alchemlyb for unified analysis
""")

# ============================================================================
# Advanced TI techniques
# ============================================================================
print("\n" + "=" * 70)
print("Advanced TI Techniques")
print("=" * 70)

print("""
1. Soft-Core Potentials:
   ---------------------
   Problem: Singularities when atoms appear/disappear

   Standard LJ: V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
                      ↑ r→0 gives infinity!

   Soft-core LJ: V_sc = 4ελ²[(σ⁶/(α(1-λ)² + r⁶))² - (σ⁶/(α(1-λ)² + r⁶))]

   Benefits:
   - Removes singularities
   - Smoother ∂H/∂λ
   - Better convergence
   - Essential for particle insertion/deletion

2. Separate Coupling:
   ------------------
   Decouple electrostatics and LJ separately:

   Stage 1: Remove charges (λ_coul: 1 → 0)
   Stage 2: Remove LJ (λ_LJ: 1 → 0)

   Or reverse for creation:
   Stage 1: Add LJ (λ_LJ: 0 → 1)
   Stage 2: Add charges (λ_coul: 0 → 1)

   Total: ΔG = ΔG_LJ + ΔG_coul

3. Dual Topology:
   --------------
   Both states A and B present, switch via λ

   V(λ) = V_common + (1-λ)V_A + λV_B

   Advantages:
   - Conceptually simple
   - Easy to implement

   Disadvantages:
   - Double the atoms
   - Need to avoid A-B interactions

4. Single Topology:
   ----------------
   Smoothly morph A → B

   Better for:
   - Small changes (mutations)
   - Same molecular topology
   - Fewer atoms

5. Hamiltonian Replica Exchange:
   -----------------------------
   Combine TI with replica exchange

   Run multiple λ values simultaneously
   Exchange between λ values

   Benefits:
   - Better sampling at each λ
   - Faster convergence
   - More robust

6. Non-equilibrium TI:
   -------------------
   Fast switching methods (Jarzynski)

   W = ∫₀¹ (∂H/∂λ) (dλ/dt) dt
   ΔG = -kT ln⟨exp(-W/kT)⟩

   Advantages:
   - Can use faster switching
   - Many short trajectories

   Disadvantages:
   - Need many realizations
   - Convergence issues
""")

# ============================================================================
# Best practices
# ============================================================================
print("\n" + "=" * 70)
print("TI Best Practices")
print("=" * 70)

print("""
1. System preparation:
   - Equilibrate state A and state B separately
   - Verify stable trajectories
   - Check energy conservation
   - Equilibrate each λ window before production

2. λ spacing:
   - Start with 11-21 windows
   - Use nonlinear spacing for soft-core
   - More points near λ=0 and λ=1
   - Test with pilot calculations
   - Refine based on ∂H/∂λ variance

3. Simulation length:
   - Minimum: 50k steps per window (small systems)
   - Typical: 200k-500k steps per window
   - Large systems: 1M+ steps per window
   - Check convergence of ⟨∂H/∂λ⟩

4. Equilibration:
   - Equilibrate 10-20% of production length
   - Start each window from previous λ
   - Monitor energy drift
   - Discard equilibration from analysis

5. Soft-core parameters:
   - Use for particle insertion/deletion
   - α = 0.5 (typical, adjustable)
   - Test different values
   - Check ∂H/∂λ smoothness

6. Error estimation:
   - Block averaging per window
   - Bootstrap analysis
   - Check time convergence
   - Compare integration methods

7. Validation:
   - Run reverse direction (hysteresis check)
   - Compare with other methods (FEP, BAR)
   - Check detailed balance
   - Validate with experimental data

8. Common pitfalls:
   - Too few λ windows → integration error
   - Poor spacing → missed features
   - Insufficient equilibration → artifacts
   - Short production → large errors
   - Wrong soft-core parameters → singularities
   - Not checking convergence → unreliable results

9. Data management:
   - Save all ∂H/∂λ trajectories
   - Keep energy files
   - Document λ schedule
   - Store analysis scripts
   - Archive configurations

10. Reporting:
    - Report all λ values used
    - Show ∂H/∂λ curve with errors
    - State integration method
    - Report simulation lengths
    - Include convergence plots
    - Compare multiple methods if possible
""")

# ============================================================================
# Practical example
# ============================================================================
print("\n" + "=" * 70)
print("Practical Example: Ligand Binding Free Energy")
print("=" * 70)

print("""
Case study: Computing absolute binding free energy

System:
  - Protein-ligand complex
  - ~30,000 atoms
  - Explicit solvent
  - Target: ΔG_bind

Thermodynamic cycle:
-------------------

     Ligand (solution) + Protein (solution)
          │                      │
          │ ΔG_bind              │
          ▼                      ▼
     Complex (solution)    Protein + Ligand (gas)
          │                      ▲
          │ ΔG_complex           │ ΔG_gas
          │                      │
          └──────────→───────────┘
              Decouple ligand

ΔG_bind = ΔG_complex - ΔG_gas

Need two TI calculations:
1. Decouple ligand in complex
2. Decouple ligand in gas phase

TI Setup:
---------

from gromos import TISimulation
import numpy as np

# Define λ schedule (separate electrostatics and LJ)
lambda_coul = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])
lambda_lj = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.0])

# Stage 1: Remove electrostatics in complex
print("Stage 1: Removing charges in complex...")
ti_coul_complex = TISimulation(
    topology_file="complex_hybrid.top",
    coordinate_file="complex_eq.cnf",
    input_file="ti_coul.imd",
    output_prefix="ti_coul_complex",
    lambda_values=lambda_coul,
    lambda_type="coulomb"
)

out_coul_complex = ti_coul_complex.run(steps_per_window=500000)
dG_coul_complex, err_coul_complex = ti_coul_complex.integrate()
print(f"ΔG_coul(complex) = {dG_coul_complex:.2f} ± {err_coul_complex:.2f} kJ/mol")

# Stage 2: Remove LJ in complex (with soft-core)
print("Stage 2: Removing LJ in complex...")
ti_lj_complex = TISimulation(
    topology_file="complex_hybrid_nocharge.top",
    coordinate_file=out_coul_complex['final_config'],
    input_file="ti_lj_softcore.imd",
    output_prefix="ti_lj_complex",
    lambda_values=lambda_lj,
    lambda_type="lj",
    soft_core=True,
    soft_core_alpha=0.5
)

out_lj_complex = ti_lj_complex.run(steps_per_window=500000)
dG_lj_complex, err_lj_complex = ti_lj_complex.integrate()
print(f"ΔG_LJ(complex) = {dG_lj_complex:.2f} ± {err_lj_complex:.2f} kJ/mol")

# Total for complex
dG_complex = dG_coul_complex + dG_lj_complex
err_complex = np.sqrt(err_coul_complex**2 + err_lj_complex**2)
print(f"ΔG_decouple(complex) = {dG_complex:.2f} ± {err_complex:.2f} kJ/mol")

# Repeat for gas phase (much faster, no solvent)
print("\\nGas phase calculations...")
dG_gas, err_gas = run_ti_gas_phase("ligand.top", "ligand.cnf")
print(f"ΔG_decouple(gas) = {dG_gas:.2f} ± {err_gas:.2f} kJ/mol")

# Final binding free energy
dG_bind = dG_complex - dG_gas
err_bind = np.sqrt(err_complex**2 + err_gas**2)

print(f"\\nFinal Result:")
print(f"ΔG_bind = {dG_bind:.2f} ± {err_bind:.2f} kJ/mol")
print(f"K_d = {np.exp(dG_bind / (0.00831446 * 300)):.2e} M")  # Assuming 300K

Expected results:
-----------------
- Total simulation time: ~50 ns × 23 windows = 1,150 ns
- Computational time: ~1-2 weeks on 128 cores
- Accuracy: ±2-5 kJ/mol (with proper sampling)
- Comparison with experiment: within error bars
""")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("""
Running TI from Python:

Quick start:
  from gromos import run_ti
  outputs = run_ti(
      topology, coordinates, input_file,
      num_windows=11,
      steps_per_window=100000,
      lambda_spacing="linear"
  )
  dG = outputs['free_energy']

Full control:
  from gromos import TISimulation
  ti = TISimulation(
      topology, coordinates, input_file,
      lambda_values=custom_lambdas
  )
  outputs = ti.run(steps_per_window=500000)
  dG, dG_err = ti.integrate(method="simpson")

Key features:
  ✓ Automatic λ scheduling
  ✓ Multiple integration methods
  ✓ Soft-core potential support
  ✓ Separate electrostatic/LJ coupling
  ✓ Comprehensive error analysis

Critical parameters:
  - Number of windows (11-51)
  - λ spacing (linear, nonlinear, optimized)
  - Steps per window (100k-1M)
  - Soft-core parameters (α ≈ 0.5)

Success criteria:
  ✓ Smooth ∂H/∂λ curve
  ✓ Low relative errors (<10% per window)
  ✓ Forward = reverse (hysteresis < 2 kJ/mol)
  ✓ Converged with simulation time

Applications:
  - Absolute binding free energies
  - Relative binding free energies
  - Solvation free energies
  - Mutation effects (ΔΔG)
  - Conformational free energies

Integration methods:
  - Trapezoid rule (simple, robust)
  - Simpson's rule (better accuracy)
  - Cubic spline (smoothest)
  - Gaussian quadrature (optimal)

Best practice workflow:
  1. Equilibrate end states
  2. Choose λ spacing (pilot runs)
  3. Equilibrate each window
  4. Run production (check convergence)
  5. Extract ∂H/∂λ
  6. Integrate (multiple methods)
  7. Error analysis
  8. Validate (reverse, comparison)

Analysis tools:
  - alchemlyb (unified interface)
  - pymbar (MBAR reweighting)
  - Custom Python scripts
  - Visualization (matplotlib)

Next steps:
  - Read TI theory papers
  - Test λ spacing on your system
  - Compare with BAR/MBAR
  - Validate with experiments
  - Consider hybrid methods (TI + BAR)
""")

print("\n" + "=" * 70)
print("Example complete!")
print("=" * 70)
print("\nKey references:")
print("  - Kirkwood, J.G. (1935) - Original TI theory")
print("  - Beutler et al. (1994) - Soft-core potentials")
print("  - Shirts & Pande (2005) - Best practices")
print("  - alchemlyb documentation - Analysis tools")
