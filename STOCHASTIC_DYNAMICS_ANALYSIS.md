# Stochastic Dynamics (Langevin) Implementation Analysis

## Files Located

### C++ Source Files
- **Header**: `/home/user/gromosXX/md++/src/algorithm/integration/stochastic.h`
- **Implementation**: `/home/user/gromosXX/md++/src/algorithm/integration/stochastic.cc`
- **Data Structure**: `/home/user/gromosXX/md++/src/topology/sd.h`
- **Random Generator**: `/home/user/gromosXX/md++/src/math/random.h`
- **Random Implementation**: `/home/user/gromosXX/md++/src/math/random.cc`

## Algorithm Overview

The GROMOS stochastic dynamics implementation uses the **Langevin equation** for integrating the equations of motion with stochastic forces. The implementation is split into four main classes:

### 1. **Stochastic_Dynamics_Vel1** (Velocity Calculation - Step 1)
- Calculates friction coefficients (gamma values) based on parameters
- Generates random vectors for stochastic forces
- Updates velocities using Langevin equations
- References equations 2.11.2.2 to 2.11.2.9 from GROMOS96 book

### 2. **Stochastic_Dynamics_Pos1** (Position Calculation - Step 1)
- Updates positions after Vel1 step
- Uses simple position update: `r_new = r_old + v_new * c6 * dt`

### 3. **Stochastic_Dynamics_Vel2** (Velocity Correction - Step 2)
- Second velocity correction (after SHAKE constraint solver if used)
- Implements equation 2.11.2.24 from GROMOS96 book
- Only runs if SHAKE constraints are enabled

### 4. **Stochastic_Dynamics_Pos2** (Position Correction - Step 2)
- Additional position correction for constrained systems
- Implements equation 2.11.2.25 and 2.11.2.26
- Uses pre-calculated random vectors from Vel1

## Core Algorithm Details

### Friction Coefficient Calculation (calc_friction_coeff)

**Input Parameters:**
- `SD`: Enable stochastic dynamics (0 or 1)
- `NTFR`: Friction coefficient mode
  - 0: Set gamma to 0.0 (no friction)
  - 1: Set gamma to CFRIC (constant friction)
  - 2: Set gamma to CFRIC × gamma₀ (per-atom scaling)
  - 3: Calculate gamma from SASA (solvent-accessible surface area)
- `CFRIC`: Global friction coefficient weighting
- `TEMP`: Bath temperature
- `NSFR`: Recalculate friction every NSFR steps
- `RCUTF`: Cutoff radius for neighbor calculation (for SASA)
- `NBREF`: Number of neighbors to consider buried (for SASA)

**Calculated Coefficients (c1-c9):**

The algorithm pre-computes 9 coefficients for each atom based on γ*dt:

#### For |γ*dt| > 0.05 (Analytical formulas):
```
gdt = gamma * dt
gdth = gdt * 0.5

emdth = exp(-gdth)
emdt = emdth²
epdt = exp(gdth)²
omdt = 1 - emdt
cdth = gdt - 3 + 4*emdth - emdt
ddth = 2 - epdth - emdth
bpdth = gdt*(epdt - 1) - 4*(epdth - 1)²
bmdth = gdt*(1 - emdt) - 4*(emdth - 1)²

c1(i) = emdt                          [exp(-gdt)]
c2(i) = omdt / gdt                   [(1 - exp(-gdt))/gdt]
c3(i) = sqrt(|omdt|)                 [sqrt(1 - exp(-gdt))]
c4(i) = sqrt(|bpdth/cdth|)           
c5(i) = gamma * ddth/cdth            
c6(i) = (epdth - emdth) / gdt        [(exp(gdt/2) - exp(-gdt/2))/gdt]
c7(i) = sqrt(|cdth|) / gamma         
c8(i) = sqrt(|bmdth/omdt|) / gamma   
c9(i) = -ddth/(gamma * omdt)         
```

#### For |γ*dt| ≤ 0.05 (Power series expansion):
Uses Taylor series expansion to avoid numerical instability when gdt is small.

### Velocity Update (Equation 2.11.2.2)

```
v_new(i) = [v_old(i) - svh] * c1(i) + F(i) * cf + vrand1(i)

where:
  cf = dt * c2(i) / m(i)
  svh = stochastic_integral_old(i) * c5(i) + vrand2(i)
  vrand1 = Gaussian random vector with σ = kT/m * c3
  vrand2 = Gaussian random vector with σ = kT/m * c4
```

### Position Update (Standard versions)

#### Vel1 → Pos1:
```
r_new(i) = r_old(i) + v_new(i) * dt * c6(i)
```

#### Vel1 → Pos2 (with constraints):
```
r_new(i) += vrand3(i) - sxh
where:
  sxh = stochastic_integral_current(i) * c9(i) + vrand4(i)
  vrand3 = Gaussian random vector with σ = kT/m * c7
  vrand4 = Gaussian random vector with σ = kT/m * c8
```

## Data Structures

### Stochastic_struct (from topology)
```cpp
struct stochastic_struct {
    math::SArray gamma;        // Friction coefficients
    math::SArray c1, c2, ..., c9;  // Pre-computed coefficients
    
    void resize(int size);
    void clear();
};
```

### RandomVectors (in Stochastic_Dynamics_Vel1)
```cpp
struct RandomVectors {
    math::VArray vrand1;  // First velocity random vector
    math::VArray vrand2;  // Second velocity random vector
    math::VArray vrand3;  // First position random vector
    math::VArray vrand4;  // Second position random vector
    
    void resize(unsigned int size);
};
```

### Configuration state tracking
```cpp
// Stochastic-specific storage
vector<Vec> stochastic_integral;     // Memory for stochastic integral
string stochastic_seed;              // RNG state for restarts
```

## Random Number Generation

### Interfaces (from random.h)
- `get()`: Uniform random in [0,1]
- `get_gauss()`: Gaussian N(0,1)
- `get_gaussian_vec()`: 3D Gaussian vector

### Implementations
1. **RandomGeneratorG96**: GROMOS96 LCG algorithm
   - Formula: `IRAND = (IRAND * 31415821 + 1) % 100000000`

2. **RandomGeneratorGSL**: GNU Scientific Library generators
   - More sophisticated algorithms
   - Can use environment variables for tuning

### Gaussian sampling
- **Box-Muller method** used internally
- Mean and standard deviation configurable

## Key Parameters and Variables

### Constants
- `k_Boltzmann = 1.987e-3 kJ/(mol·K)` or `8.314462618e-3 kJ/(mol·K)`

### Per-atom storage requirements
- gamma: 1 double
- c1-c9: 9 doubles
- stochastic_integral: 1 Vec3 (3 floats)
- **Total**: ~13 KB per 1000 atoms

## Boundary Conditions

The implementation uses templated boundary condition handling:
- **Vacuum**: No periodic boundary conditions
- **Rectangular**: Rectangular periodic box
- **Truncated Octahedron**: Non-orthogonal cell
- Used for neighbor calculation in SASA mode (NTFR=3)

## Integration Sequence in GROMOS

A typical MD simulation cycle:
1. `StochasticDynamics_Vel1::apply()` → Calculate velocities, generate random vectors
2. `StochasticDynamics_Pos1::apply()` → Update positions
3. Force calculation
4. (Optional) `SHAKE` constraint solver
5. (Optional) `StochasticDynamics_Vel2::apply()` → Velocity correction
6. (Optional) `StochasticDynamics_Pos2::apply()` → Position correction

## Important Notes

### Energy Conservation
- Uses specific coefficients (c1-c9) derived from exact solutions to Langevin equation
- Two algorithm paths: analytical (|gdt| > 0.05) vs power series (|gdt| ≤ 0.05)

### Friction Modes
- **SASA-based** (NTFR=3): Calculates neighbor count within RCUTF
  - Useful for protein simulations (bury heavy atoms less)
  - Formula: `gamma = cfric * max(0, 1 - neigh/nbref)`

### Constraints Compatibility
- Supports SHAKE constraints (LINCS not mentioned as compatible)
- Vel2/Pos2 steps correct velocities/positions after constraint application

### Initial Stochastic Integral
- Can be generated automatically or read from trajectory file
- Uses same RNG seed for reproducibility
- Stored in configuration state

### Parallel Efficiency
- Friction calculation is O(N²) for NTFR=3 (neighbor counting)
- Should recalculate only periodically (NSFR)
- Velocity/position updates are O(N) and easily parallelizable

## References in Code
- References to GROMOS96 manual equations (2.11.2.x)
- Based on Langevin dynamics from van Gunsteren & Berendsen

---

## Rust Port Considerations

### Key Components to Port
1. Friction coefficient pre-computation (analytical and series paths)
2. Gaussian random vector generation
3. Velocity update equations (with numerical stability checks)
4. Position update equations
5. Neighbor-based SASA friction calculation
6. State persistence (stochastic_integral, seed)

### Challenge Areas
- **Numerical precision**: Power series expansion for small gdt values
- **Random number reproducibility**: Must maintain seed compatibility
- **Memory layout**: Pre-computed coefficients stored per-atom

### SIMD Opportunities
- Parallel friction coefficient calculation for independent atoms
- Batch Gaussian random number generation
- Parallel velocity/position updates

