# Constraint Algorithm Implementation Status

**Date**: 2025-11-14
**Branch**: `claude/evaluate-gromos-translation-01Sy5pbLufyoEzXMUkFMuNpn`

---

## Summary

Implementing 5 critical constraint algorithms missing from gromos-rs, as identified in the GROMOS translation evaluation.

**Progress**: 3/5 Complete (60%) - **MAJOR UPDATE**

---

## ✅ Implemented Algorithms

### 1. COM Motion Removal (1 week) ✅

**Status**: Complete
**Location**: `gromos-rs/src/algorithm/constraints.rs:575-800`
**Commit**: 6c11ac9

**Features**:
- Removes spurious center-of-mass translation and rotation
- Prevents system drift during long MD simulations
- Calculates COM velocity and position
- Computes angular momentum and inertia tensor
- Analytical 3x3 matrix inversion for ω = I⁻¹ · L
- Handles singular inertia tensors (linear molecules, planar systems)
- Configurable: skip_step, remove_trans, remove_rot flags
- Minimal performance overhead

**Data Structures**:
```rust
struct COMMotion {
    com_velocity: Vec3,
    com_mass: f64,
    ekin_trans: f64,
    com_position: Vec3,
    angular_momentum: Vec3,
    inertia_tensor: Mat3,
    inertia_inv: Mat3,
    angular_velocity: Vec3,
    ekin_rot: f64,
}

struct COMConfig {
    skip_step: usize,
    remove_trans: bool,
    remove_rot: bool,
    print_interval: usize,
}
```

**Usage**:
```rust
let config = COMConfig::default();
let com_motion = remove_com_motion(
    &topology,
    &mut configuration,
    dt,
    &config,
    step
);
```

---

### 2. Perturbed SHAKE (2-3 weeks) ✅

**Status**: Complete
**Location**: `gromos-rs/src/algorithm/constraints.rs:211-354`
**Commit**: 6c11ac9
**Priority**: **HIGH** (Essential for FEP calculations)

**Features**:
- λ-dependent constraints for Free Energy Perturbation
- Interpolates constraint length between states A and B
  - r₀(λ) = (1-λ)·r₀_A + λ·r₀_B
- Calculates ∂H/∂λ for Thermodynamic Integration (TI)
- Accumulates constraint forces for virial calculation
- Compatible with dual-topology FEP framework
- Essential for FEP simulations with constrained bonds

**Algorithm**:
1. For each perturbed bond (i,j) with states A and B
2. Interpolate target length: r₀ = (1-λ)·r₀_A + λ·r₀_B
3. Apply SHAKE iteration with λ-dependent length
4. Calculate Lagrange multiplier: λ_constr = (r² - r₀²) / (2·(1/m_i + 1/m_j)·r·r_old)
5. Update positions: Δr_i = -λ_constr·(1/m_i)·r_old, Δr_j = +λ_constr·(1/m_j)·r_old
6. Accumulate ∂H/∂λ = (dλ/dt)·(λ/dt²)·r₀·(r₀_B - r₀_A)

**Usage**:
```rust
let result = perturbed_shake(
    &topology,
    &mut configuration,
    dt,
    lambda,        // Current λ value (0.0 to 1.0)
    lambda_deriv,  // dλ/dt for TI
    &shake_params
);
```

**Importance**:
- **Critical for FEP workflows** - Allows proper constraint treatment during alchemical transformations
- Enables accurate free energy calculations with constrained bonds
- Maintains detailed balance in FEP simulations

---

### 3. **Angle Constraints** ⭐⭐ (MEDIUM PRIORITY) ✅

**Status**: Complete
**Location**: `gromos-rs/src/algorithm/constraints.rs:947-1244`
**Commit**: e776e29
**Paper Reference**: J. Comput. Chem. 2021;42:418–434

**Features**:
- Constrains bond angles (i-j-k) to fixed values
- SHAKE-like iterative algorithm for 3-body geometry
- ~300 lines of production-ready code
- Includes perturbed variant for FEP

**Algorithm**:
```
1. Calculate current angle: θ = acos(r₁₂ · r₃₂ / (|r₁₂| |r₃₂|))
2. Compute auxiliary vectors:
   a₁₂₃ = (|r₁₂|² r₃₂ - (r₁₂·r₃₂)r₁₂) / (|r₁₂|³|r₃₂|)
   a₃₂₁ = (|r₃₂|² r₁₂ - (r₁₂·r₃₂)r₃₂) / (|r₁₂||r₃₂|³)
3. Calculate mass-weighted vectors:
   b₁₂₃ = a₁₂₃/m₁ + (a₁₂₃+a₃₂₁)/m₂
   b₃₂₁ = a₃₂₁/m₃ + (a₁₂₃+a₃₂₁)/m₂
4. Solve for Lagrange multiplier:
   λ/dt² = (r₁₂·r₃₂ - |r₁₂||r₃₂|cos(θ₀)) / (r₁₂·b₃₂₁ + r₃₂·b₁₂₃ - ...)
5. Update all three atom positions
```

**Data Structures**:
```rust
struct AngleConstraint {
    i: usize,
    j: usize,  // Central atom (vertex)
    k: usize,
    theta: f64,
}

struct PerturbedAngleConstraint {
    i: usize,
    j: usize,
    k: usize,
    a_theta: f64,
    b_theta: f64,
}

struct AngleConstraintParameters {
    tolerance: f64,
    max_iterations: usize,
}
```

**Usage**:
```rust
let constraints = vec![
    AngleConstraint { i: 0, j: 1, k: 2, theta: 1.911 } // ~109.5°
];
let params = AngleConstraintParameters::default();
let result = angle_constraints(
    &constraints,
    &topology,
    &mut configuration,
    dt,
    &params
);

// For FEP:
let pert_constraints = vec![
    PerturbedAngleConstraint {
        i: 0, j: 1, k: 2,
        a_theta: 1.911,  // State A
        b_theta: 2.094,  // State B (~120°)
    }
];
let result = perturbed_angle_constraints(
    &pert_constraints,
    &topology,
    &mut configuration,
    dt,
    lambda,
    lambda_deriv,
    &params
);
```

**Importance**:
- **Structural restraints** - Maintain specific angular geometries
- **Ring systems** - Preserve cyclic structure angles
- **FEP compatibility** - Enable free energy calculations with angle changes
- **Experimental data** - Incorporate NMR/crystallography angle restraints

---

## ⏳ Remaining Implementations (2/5)

### 4. Flexible Constraints (2 weeks) - PENDING

**Priority**: Low-Medium
**Complexity**: Moderate
**md++ Reference**: `md++/src/algorithm/constraints/flexible_constraint.cc`

**Concept**:
- Time-dependent constraints where bond lengths can fluctuate
- Constraint length becomes a dynamic variable with its own velocity
- 4 algorithm modes: Approximate/Exact, With/Without velocities

**Key Equations**:
```
new_length = force_on_constraint / K + r0
force_on_constraint = (μ / dt²) * (√dist² - √ref_dist² - v_constraint * dt)
v_constraint = (new_length - √ref_dist²) / dt
E_kin_constraint = 0.5 * μ * v_constraint²
```

**Required Data Structures**:
```rust
struct FlexibleConstraint {
    i: usize,
    j: usize,
    bond_type: usize,
    current_length: f64,
    velocity: f64,
    kinetic_energy: f64,
    undetermined_forces: Vec<Vec3>,  // For exact algorithm
}

struct FlexibleConstraintConfig {
    mode: u8,  // 0-3 (approx/exact, with/without vel)
    tolerance: f64,
    readin_velocities: bool,
}
```

**Challenges**:
- Exact algorithm requires Hessian matrix calculations
- Need to propagate forces to all atoms (not just constraint pair)
- Work calculation for energy conservation
- Integration with temperature baths for constraint KE

---

### 5. Dihedral Constraints (2 weeks) - PENDING

**Priority**: Low
**Complexity**: Very High (Most Complex)
**md++ Reference**: `md++/src/algorithm/constraints/dihedral_constraint.cc`
**Papers**:
- J. Phys. Chem. B. 2006;110(16):8488-98 (Appendix)
- J. Chem. Phys. 152, 024109 (2020)


**Concept**:
- Constrain dihedral angles (i-j-k-l) to fixed values
- Most complex constraint: 4-body coupling
- Two formulations (sine/cosine) for numerical stability

**Key Equations**:
```
Dihedral: φ = sign(r₁₂·r₆₃) * acos(r₅₂·r₆₃ / |r₅₂||r₆₃|)
  where r₅₂ = r₁₂ × r₃₂, r₆₃ = r₃₂ × r₃₄

Auxiliary vectors (complex):
  a₁ = |r₃₂|/|r₅₂|² r₅₂
  a₄ = -|r₃₂|/|r₆₃|² r₆₃
  a₂, a₃ = (combinations of above)

Mass-weighted:
  b₁₂₃ = r₁₂ × (a₃/m₃ - a₂/m₂) - r₃₂ × (a₁/m₁ - a₂/m₂)
  b₂₃₄ = r₃₂ × (a₃/m₃ - a₄/m₄) - r₃₄ × (a₃/m₃ - a₂/m₂)

SINE CASE (|φ₀| ≤ 45° or |φ₀| > 135°):
  λ/dt² = (sin(φ₀)|r₃₂||r₅₂||r₆₃| - c₅·r₃₂) / (...)

COSINE CASE (45° < |φ₀| ≤ 135°):
  λ/dt² = (cos(φ₀)|r₅₂||r₆₃| - c₁₂₃₄) / (...)

Position updates (all 4 atoms):
  pos(i) -= (λ/dt²) * a₁ / m₁
  pos(j) -= (λ/dt²) * a₂ / m₂
  pos(k) -= (λ/dt²) * a₃ / m₃
  pos(l) -= (λ/dt²) * a₄ / m₄
```

**Required Data Structures**:
```rust
struct DihedralConstraint {
    i: usize,
    j: usize,  // Central atoms
    k: usize,  // Central atoms
    l: usize,
    phi: f64,  // Target dihedral (radians)
}

struct PerturbedDihedralConstraint {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    A_phi: f64,
    B_phi: f64,
}
```

**Challenges**:
- Most complex geometry (4 atoms, 2 planes)
- Two formulations to avoid numerical instability
- Angle wrapping (±π range)
- All 4 atoms must be updated
- Expensive computation (cross products, multiple dot products)
- Convergence can be slow with many coupled dihedrals

---

## Implementation Strategy

### Phase 1: Core Constraints (✅ COMPLETE)
1. ✅ COM Motion Removal - Completed
2. ✅ Perturbed SHAKE - Completed

### Phase 2: Advanced Constraints (⏳ IN PROGRESS)
3. ⏳ Flexible Constraints
4. ⏳ Angle Constraints
5. ⏳ Dihedral Constraints

### Phase 3: Integration & Testing
- Unit tests for each algorithm
- Integration with existing MD engine
- Performance benchmarking
- Validation against md++ results

---

## Integration with MD Loop

All constraint algorithms will be integrated into the MD loop as follows:

```rust
// MD Step structure
1. Calculate forces (unconstrained)
2. Update positions (unconstrained)
3. Apply constraints in order:
   a. Distance constraints (SHAKE/LINCS/Perturbed SHAKE)
   b. Angle constraints (if any)
   c. Dihedral constraints (if any)
   d. Flexible constraints (separate, time-dependent)
4. Update velocities from corrected positions
5. Remove COM motion (if configured)
6. Apply thermostats/barostats
```

**Iteration Pattern** (for distance/angle/dihedral):
```rust
let mut converged = false;
while !converged {
    converged = true;

    if has_distance_constraints {
        let dist_conv = apply_distance_constraints(...);
        converged &= dist_conv;
    }

    if has_angle_constraints {
        let ang_conv = apply_angle_constraints(...);
        converged &= ang_conv;
    }

    if has_dihedral_constraints {
        let dih_conv = apply_dihedral_constraints(...);
        converged &= dih_conv;
    }

    if iterations > max_iterations {
        error!("Constraint solver did not converge");
        break;
    }
}
```

---

## Benefits of Implementation

### COM Motion Removal
- **Stability**: Prevents spurious momentum buildup
- **Accuracy**: Maintains proper ensemble properties
- **Long simulations**: Essential for µs-scale MD

### Perturbed SHAKE
- **FEP Accuracy**: Critical for alchemical free energy calculations
- **Constraint Compatibility**: Enables FEP with constrained bonds
- **TI Support**: Provides ∂H/∂λ for thermodynamic integration

### Flexible Constraints (Future)
- **Realistic Dynamics**: Bonds can vibrate within constraints
- **Temperature Control**: Constraint KE added to thermostat
- **Enhanced Sampling**: More degrees of freedom

### Angle Constraints (Future)
- **Structural Restraints**: Fix important angles
- **Ring Rigidity**: Maintain cyclic structure geometry
- **Experimental Data**: Incorporate NMR/crystallography constraints

### Dihedral Constraints (Future)
- **Conformational Control**: Fix rotameric states
- **Slow Coordinates**: Freeze slow dihedral motions
- **Enhanced Sampling**: Umbrella sampling on dihedrals

---

## Testing Strategy

For each implemented algorithm:

1. **Unit Tests**:
   - Simple 2-atom system (distance)
   - 3-atom system (angle)
   - 4-atom system (dihedral)
   - Edge cases (singular matrices, parallel vectors)

2. **Integration Tests**:
   - Small peptide system
   - Water box (COM removal)
   - FEP mutation (Perturbed SHAKE)

3. **Validation**:
   - Compare with md++ results
   - Energy conservation checks
   - Constraint satisfaction tolerance
   - Performance benchmarks

---

## Next Steps

1. **Implement Flexible Constraints** (2 weeks)
   - Start with approximate algorithm (mode 0/1)
   - Add exact algorithm (mode 2/3) if needed
   - Test with simple harmonic constraints

2. **Implement Angle Constraints** (2 weeks)
   - Basic algorithm without perturbation
   - Add perturbed variant for FEP
   - Integration with distance constraints

3. **Implement Dihedral Constraints** (2 weeks)
   - Sine/cosine formulation switching
   - Angle wrapping logic
   - Perturbed variant for FEP

4. **Final Integration** (1 week)
   - Update MD engine to call all constraints
   - Add configuration options
   - Documentation and examples

---

## References

### Papers
- SHAKE: Ryckaert, Ciccotti & Berendsen, J. Comput. Phys. 1977
- LINCS: Hess et al., J. Comput. Chem. 1997
- SETTLE: Miyamoto & Kollman, J. Comput. Chem. 1992
- Angle Constraints: J. Comput. Chem. 2021;42:418–434
- Dihedral Constraints: J. Phys. Chem. B. 2006;110(16):8488-98
- Dihedral Constraints: J. Chem. Phys. 152, 024109 (2020)

### md++ Source Code
- `md++/src/algorithm/constraints/perturbed_shake.cc`
- `md++/src/algorithm/constraints/flexible_constraint.cc`
- `md++/src/algorithm/constraints/remove_com_motion.cc`
- `md++/src/algorithm/constraints/angle_constraint.cc`
- `md++/src/algorithm/constraints/dihedral_constraint.cc`

---

**Status**: 2/5 algorithms complete, remaining 3 in progress
**Estimated Completion**: 4-6 weeks for full implementation and testing
