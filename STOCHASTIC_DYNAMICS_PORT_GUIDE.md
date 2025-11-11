# Stochastic Dynamics Langevin Implementation - Rust Port Guide

## Quick Reference

This guide consolidates the stochastic dynamics (Langevin) implementation from GROMOS C++ for porting to Rust.

## Files Created in This Analysis

1. **STOCHASTIC_DYNAMICS_ANALYSIS.md** - Complete algorithm overview with parameter descriptions
2. **STOCHASTIC_DYNAMICS_CODE_REFERENCE.md** - Key C++ code snippets with line numbers
3. **STOCHASTIC_DYNAMICS_PORT_GUIDE.md** - This file

## Original C++ Source Files

### Main Implementation
- `/home/user/gromosXX/md++/src/algorithm/integration/stochastic.h` - Class definitions
- `/home/user/gromosXX/md++/src/algorithm/integration/stochastic.cc` - Implementation (~528 lines)

### Supporting Files
- `/home/user/gromosXX/md++/src/topology/sd.h` - Stochastic data structure
- `/home/user/gromosXX/md++/src/math/random.h` - Random number generator interface
- `/home/user/gromosXX/md++/src/math/random.cc` - RNG implementations
- `/home/user/gromosXX/md++/src/simulation/parameter.h` (lines 2980-3035) - Parameter definitions

## Algorithm Summary

### Four Main Classes

1. **Stochastic_Dynamics_Vel1**
   - Calculates friction coefficients (gamma)
   - Generates 4 random vectors per atom
   - Updates velocities (Langevin dynamics)
   - Location: stochastic.cc:207-423

2. **Stochastic_Dynamics_Pos1**
   - Updates positions after velocity calculation
   - Location: stochastic.cc:429-446

3. **Stochastic_Dynamics_Vel2**
   - Corrects velocities for SHAKE constraints
   - Location: stochastic.cc:452-499

4. **Stochastic_Dynamics_Pos2**
   - Corrects positions for SHAKE constraints
   - Location: stochastic.cc:505-527

### Critical Algorithms

#### 1. Friction Coefficient Calculation (calc_friction_coeff)
- Location: stochastic.cc:49-202
- Two paths: Analytical (|gdt| > 0.05) and Power Series (|gdt| ≤ 0.05)
- Calculates 9 coefficients (c1-c9) per atom
- Supports 4 friction modes:
  - NTFR=0: gamma = 0 (no friction)
  - NTFR=1: gamma = CFRIC (constant)
  - NTFR=2: gamma = CFRIC × gamma₀ (per-atom)
  - NTFR=3: gamma from SASA (neighbor-based)

#### 2. Velocity Update (Equation 2.11.2.2)
- Location: stochastic.cc:333-415
- Uses pre-computed coefficients
- Applies random forces from Maxwell-Boltzmann distribution
- Formula: `v_new = (v_old - svh) * c1 + F/m * cf + vrand1`

#### 3. Position Update (Pos1)
- Location: stochastic.cc:439-440
- Simple: `r_new = r_old + v_new * dt * c6`

#### 4. Constraint Corrections (Vel2/Pos2)
- Location: stochastic.cc:481-487 (Vel2), 517-521 (Pos2)
- Only applied if SHAKE is active

## Key Parameters

### Input Parameters (from simulation.parameter)
```
SD          Enable stochastic dynamics (0/1)
NTFR        Friction coefficient mode (0-3)
CFRIC       Global friction weighting (1/ps)
TEMP        Bath temperature (K)
NSFR        Recalculate friction every N steps
RCUTF       Neighbor cutoff (nm)
NBREF       Buried neighbor threshold (count)
```

### Per-Atom State
```
gamma              Friction coefficient (1/ps)
c1 - c9            Pre-computed coefficients
stochastic_integral Vec3 memory for damping term
```

## Numerical Requirements

### Constants
- Boltzmann constant: 8.314462618 kJ/(mol·K)

### Stability Thresholds
- Analytical formula used when: |gamma * dt| > 0.05
- Power series used when: |gamma * dt| ≤ 0.05
- Defensive sqrts: Always use abs() before sqrt()

### Precision
- Double precision required throughout
- Float acceptable for velocity/position storage
- Need stable Gaussian random number generation

## Memory Requirements

Per atom:
- gamma: 8 bytes
- c1-c9: 72 bytes (9 × 8 bytes)
- stochastic_integral: 12 bytes (3 × 4 bytes)
- **Total**: ~92 bytes per atom

For 1000 atoms: ~90 KB
For 100,000 atoms: ~9 MB

## Rust Implementation Checklist

- [ ] Create `src/algorithm/stochastic_dynamics.rs` module
- [ ] Implement `calc_coefficients_analytical()`
- [ ] Implement `calc_coefficients_power_series()`
- [ ] Implement SASA neighbor calculation
- [ ] Implement Gaussian RNG (using rand crate or custom)
- [ ] Implement velocity update loop
- [ ] Implement position update loop
- [ ] Add Vel2/Pos2 constraint correction paths
- [ ] Create comprehensive tests comparing against C++
- [ ] Verify energy conservation in NVE simulations
- [ ] Test reproducibility with fixed seed
- [ ] Validate numerical accuracy

## Integration with Existing Rust Code

### Proposed Trait Implementation
```rust
impl Integrator for StochasticDynamics {
    fn step(&mut self, dt: f64, topo: &Topology, conf: &mut Configuration) {
        // Vel1 + Pos1 or Vel2 + Pos2 depending on constraints
    }
}
```

### Data Structure (SoA recommended)
```rust
pub struct StochasticState {
    pub gamma: Vec<f64>,
    pub c1: Vec<f64>, c2: Vec<f64>, ..., c9: Vec<f64>,
    pub stochastic_integral: Vec<Vec3>,
    pub rng_seed: String,
}
```

## Testing Strategy

1. **Unit tests**: Test each coefficient calculation function
2. **Reference tests**: Compare with C++ output on small systems
3. **Conservation tests**: Verify energy conservation in NVE
4. **Reproducibility tests**: Same seed → same trajectory
5. **Edge cases**: Test with gdt near 0.05 threshold
6. **SASA tests**: Verify neighbor counting

## References

- GROMOS96 manual (equations 2.11.2.x)
- Schmid et al., "Architecture, Implementation and Parallelisation of GROMOS", 
  Comp. Phys. Commun. 183 (2012)
- van Gunsteren & Berendsen, Langevin dynamics theory

## Performance Notes

- Friction calculation O(N²) for NTFR=3 (should recalculate infrequently)
- Velocity/position updates are O(N)
- 4 Gaussian random vectors generated per atom per step
- Excellent parallelization potential with Rayon

## Known Issues in C++ (for Rust avoidance)

1. **Neighbor calculation**: O(N²) loop - consider spatial decomposition for large systems
2. **Random number generation**: GSL vs G96 differences possible
3. **Numerical precision**: Power series can have accumulated errors
4. **Constraint interaction**: Only tested with SHAKE, not LINCS

## Implementation Priority

1. **High**: Velocity update equations, coefficient calculation
2. **Medium**: SASA friction calculation, constraint corrections
3. **Low**: Alternative RNG implementations (use rand crate for now)

---

**Last Updated**: 2025-11-10
**Status**: Analysis complete, ready for implementation
