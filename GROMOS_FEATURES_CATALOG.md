# GROMOS md++ Complete Features Catalog

**Date**: 2025-11-10
**Source**: Analyzed 486 C++ files in md++/src
**Purpose**: Comprehensive inventory of all GROMOS simulation capabilities

---

## Executive Summary

GROMOS is **far more than a basic MD engine**. It includes:
- **13 integration algorithms** (not just leap-frog)
- **23+ special interactions** (restraints, QM/MM, X-ray)
- **9 constraint algorithms** (SHAKE, SETTLE, LINCS, etc.)
- **Enhanced sampling** methods (replica exchange, metadynamics, EDS)
- **Free energy** calculations (TI, slow growth, perturbation)
- **Structure refinement** with experimental data (NMR, X-ray)
- **Advanced analysis** capabilities

---

## 1. Integration Algorithms (13 types)

### 1.1 Standard Molecular Dynamics
| Algorithm | File | Purpose | Rust Status |
|-----------|------|---------|-------------|
| **Leap-Frog** | `leap_frog.cc` | Standard MD integrator | ✅ Implemented |
| **Scaled Leap-Frog** | `scaled_leap_frog.cc` | Multiple time-stepping | ❌ Not implemented |
| **Stochastic Dynamics** | `stochastic.cc` | Langevin/BD dynamics | ❌ Not implemented |

### 1.2 Energy Minimization
| Algorithm | File | Purpose | Rust Status |
|-----------|------|---------|-------------|
| **Steepest Descent** | `steepest_descent.cc` | Simple minimization | ❌ Not implemented |
| **Conjugate Gradient** | `conjugate_gradient.cc` | Efficient minimization | ❌ Not implemented |

### 1.3 Enhanced Sampling
| Algorithm | File | Purpose | Rust Status |
|-----------|------|---------|-------------|
| **Monte Carlo** | `monte_carlo.cc` | Chemical MC sampling | ❌ Not implemented |
| **EDS** | `eds.cc` | Enveloping Distribution Sampling | ❌ Not implemented |
| **GaMD** | `gamd.cc` | Gaussian Accelerated MD | ❌ Not implemented |
| **Multigradient** | `multigradient.cc` | Multiple potential sampling | ❌ Not implemented |

### 1.4 Free Energy Methods
| Algorithm | File | Purpose | Rust Status |
|-----------|------|---------|-------------|
| **Slow Growth** | `slow_growth.cc` | TI with continuous λ | ❌ Not implemented |
| **Lattice Shift** | `lattice_shift.cc` | Lattice sum derivatives | ❌ Not implemented |

### 1.5 Analysis & Special
| Algorithm | File | Purpose | Rust Status |
|-----------|------|---------|-------------|
| **Energy Calculation** | `energy_calculation.cc` | Single-point energies | ✅ Partially (tests) |
| **Analyze** | `analyze.cc` | Trajectory analysis | ❌ Not implemented |

**Key Insight**: gromos-rs currently only has **leap-frog**. Need 12 more integration methods!

---

## 2. Constraint Algorithms (9 types)

### 2.1 Bond Constraints
| Algorithm | File | Purpose | Implementation |
|-----------|------|---------|----------------|
| **SHAKE** | `shake.cc` | Iterative bond constraints | Standard algorithm |
| **M-SHAKE** | `m_shake.cc` | Mass-weighted SHAKE | Matrix formulation |
| **Perturbed SHAKE** | `perturbed_shake.cc` | SHAKE for λ-perturbation | Free energy |
| **GPU SHAKE** | `gpu_shake.cc` | CUDA-accelerated SHAKE | GPU version |
| **LINCS** | `lincs.cc` | Linear constraint solver | Alternative to SHAKE |

### 2.2 Special Constraints
| Algorithm | File | Purpose | Implementation |
|-----------|------|---------|----------------|
| **SETTLE** | `settle.cc` | Analytical water constraints | Fast rigid water |
| **GPU SETTLE** | `gpu_settle.cc` | CUDA SETTLE | GPU version |
| **Flexible Constraints** | `flexible_constraint.cc` | Time-dependent constraints | Dynamic constraints |
| **Position Constraints** | `position_constraints.cc` | Fix atom positions | Harmonic restraints |

### 2.3 System Constraints
| Algorithm | File | Purpose | Implementation |
|-----------|------|---------|----------------|
| **COM Motion Removal** | `remove_com_motion.cc` | Remove drift | System cleanup |
| **Rotation/Translation** | `rottrans.cc` | Remove rotation | System cleanup |
| **Angle Constraints** | `angle_constraint.cc` | Fix angles | Rare, special cases |
| **Dihedral Constraints** | `dihedral_constraint.cc` | Fix dihedrals | Rare, special cases |

**Rust Status**: ❌ **No constraint algorithms implemented yet**

---

## 3. Bonded Interactions (14 types)

### 3.1 Standard Force Field
| Interaction | File | Formula | Status |
|-------------|------|---------|--------|
| **Quartic Bonds** | `quartic_bond_interaction.cc` | GROMOS quartic | ⚠️ Energy only |
| **Harmonic Bonds** | `harmonic_bond_interaction.cc` | Simple harmonic | ❌ Not implemented |
| **CG Bonds** | `cg_bond_interaction.cc` | Coarse-grained | ❌ Not implemented |
| **Angles** | `angle_interaction.cc` | Cosine-based | ⚠️ Energy only |
| **Harmonic Angles** | `harm_angle_interaction.cc` | Simple harmonic | ❌ Not implemented |
| **Proper Dihedrals** | `dihedral_interaction.cc` | Torsional | ❌ Not implemented |
| **New Dihedrals** | `dihedral_new_interaction.cc` | Improved formula | ❌ Not implemented |
| **Improper Dihedrals** | `improper_dihedral_interaction.cc` | Out-of-plane | ❌ Not implemented |
| **Cross-Dihedrals** | `crossdihedral_interaction.cc` | 8-atom term | ❌ Not implemented |

### 3.2 Perturbed (Free Energy)
All bonded terms have perturbed versions for λ-dependent free energy calculations:
- `perturbed_quartic_bond_interaction.cc`
- `perturbed_harmonic_bond_interaction.cc`
- `perturbed_angle_interaction.cc`
- `perturbed_dihedral_interaction.cc`
- `perturbed_improper_dihedral_interaction.cc`

Plus soft-core variants for avoiding singularities:
- `perturbed_soft_bond_interaction.cc`
- `perturbed_soft_angle_interaction.cc`
- `perturbed_soft_improper_interaction.cc`

**Status**: ❌ **No perturbed interactions in gromos-rs**

---

## 4. Nonbonded Interactions

### 4.1 Basic Nonbonded
- **Lennard-Jones** - ✅ Implemented in gromos-rs
- **Coulomb electrostatics** - ⚠️ Partially implemented
- **Reaction field** - ⚠️ Partially implemented

### 4.2 Long-Range Electrostatics
- **Ewald summation**
- **Particle Mesh Ewald (PME)**
- **P3M (Particle-Particle Particle-Mesh)**
- **Lattice sum**

**Status**: ❌ **No long-range electrostatics in gromos-rs**

### 4.3 Pairlist Algorithms
| Algorithm | File | Method | Status |
|-----------|------|--------|--------|
| **Standard** | `standard_pairlist_algorithm.cc` | Simple distance | ⚠️ Basic version |
| **Grid Cell** | `grid_cell_pairlist.cc` | Spatial decomposition | ❌ Not implemented |
| **Extended Grid** | `extended_grid_pairlist_algorithm.cc` | Optimized grid | ❌ Not implemented |

---

## 5. Special Interactions (23+ types)

### 5.1 Experimental Restraints (NMR)
| Restraint | File | Experimental Data | Purpose |
|-----------|------|-------------------|---------|
| **Distance** | `distance_restraint_interaction.cc` | NOE | Structure refinement |
| **Angle** | `angle_restraint_interaction.cc` | Angular constraints | Backbone angles |
| **Dihedral** | `dihedral_restraint_interaction.cc` | Torsion angles | φ/ψ angles |
| **J-value** | `jvalue_restraint_interaction.cc` | NMR J-coupling | Through-bond coupling |
| **RDC** | `rdc_restraint_interaction.cc` | Residual dipolar coupling | Orientation |
| **Order Parameter** | `order_parameter_restraint_interaction.cc` | S² order parameters | Dynamics |

**Purpose**: Refine MD structures using experimental NMR data

**Status**: ❌ **None implemented in gromos-rs**

### 5.2 X-ray Crystallography
| Feature | File | Purpose | Status |
|---------|------|---------|--------|
| **X-ray Restraints** | `xray_restraint_interaction.cc` | Structure factors | ❌ Not implemented |
| **X-ray Density** | `special/xray/dens.cc` | Electron density | ❌ Not implemented |
| **Structure Factors** | `special/xray/sf.cc` | Diffraction calculation | ❌ Not implemented |

**Purpose**: Refine MD structures against X-ray diffraction data

### 5.3 Enhanced Sampling
| Method | File | Purpose | Status |
|--------|------|---------|--------|
| **Local Elevation** | `local_elevation_interaction.cc` | Metadynamics-like | ❌ Not implemented |
| **Distance Field** | `distance_field_interaction.cc` | Biasing potential | ❌ Not implemented |
| **B&S Potential** | `bs_interaction.cc` | Boundary smoothing | ❌ Not implemented |
| **ADDE Reweighting** | `adde_reweighting.cc` | Adaptive sampling | ❌ Not implemented |

### 5.4 External Forces
| Feature | File | Purpose | Status |
|---------|------|---------|--------|
| **Electric Field** | `electric_field_interaction.cc` | Applied E-field | ❌ Not implemented |
| **NEMD** | `nemd.cc` | Non-equilibrium MD | ❌ Not implemented |
| **Pressure Scaling** | `pscale.cc` | Anisotropic pressure | ❌ Not implemented |

### 5.5 Symmetry & Positioning
| Feature | File | Purpose | Status |
|---------|------|---------|--------|
| **Position Restraints** | `position_restraint_interaction.cc` | Fix/restrain atoms | ❌ Not implemented |
| **Symmetry Restraints** | `symmetry_restraint_interaction.cc` | Enforce symmetry | ❌ Not implemented |

### 5.6 Free Energy Variants
All restraints have perturbed (λ-dependent) versions:
- `perturbed_distance_restraint_interaction.cc`
- `perturbed_angle_restraint_interaction.cc`
- `perturbed_dihedral_restraint_interaction.cc`
- `perturbed_distance_field_interaction.cc`
- `eds_distance_restraint_interaction.cc` (EDS-specific)

---

## 6. QM/MM Hybrid Methods

### 6.1 Quantum Mechanics Integration
| Feature | File | Purpose | Status |
|---------|------|---------|--------|
| **QM/MM Interaction** | `qmmm/qmmm_interaction.cc` | Hybrid QM/MM | ❌ Not implemented |
| **QM Worker** | `qmmm/qm_worker.cc` | External QM interface | ❌ Not implemented |
| **NN Worker** | `qmmm/nn_worker.cc` | Neural network QM | ❌ Not implemented |
| **QM/MM Nonbonded** | `qmmm/nonbonded/` | QM-MM interactions | ❌ Not implemented |

**Purpose**: Treat critical region (e.g., active site) with quantum mechanics, rest with MM

**Supported QM Programs**:
- Turbomole
- MNDO
- Gaussian
- MOPAC
- Neural network potentials

**Status**: ❌ **Completely missing in gromos-rs** (advanced feature)

---

## 7. Thermostats & Barostats

### 7.1 Temperature Control
| Method | File | Type | Status |
|--------|------|------|--------|
| **Berendsen** | `temperature/thermostat.cc` | Weak coupling | ❌ Not implemented |
| **Nosé-Hoover** | `temperature/thermostat.cc` | Extended system | ❌ Not implemented |
| **Andersen** | `stochastic.cc` | Stochastic collisions | ❌ Not implemented |

### 7.2 Pressure Control
| Method | File | Type | Status |
|--------|------|------|--------|
| **Berendsen** | `pressure/berendsen_barostat.cc` | Weak coupling | ❌ Not implemented |
| **Parrinello-Rahman** | `pressure/` | Full anisotropic | ❌ Not implemented |

**Status**: ❌ **No thermostats/barostats in gromos-rs**

---

## 8. Replica Exchange Methods

### 8.1 Standard REMD
| Type | Directory | Purpose | Status |
|------|-----------|---------|--------|
| **Temperature REMD** | `replicaExchange/` | T-REMD | ❌ Not implemented |
| **Hamiltonian REMD** | `replicaExchange/` | H-REMD | ❌ Not implemented |
| **2D REMD** | `replica_exchangers/2D_T_lambda_REPEX/` | T+λ exchange | ❌ Not implemented |

### 8.2 EDS Replica Exchange
| Type | Directory | Purpose | Status |
|------|-----------|---------|--------|
| **1D EDS-RE** | `replica_exchangers/1D_S_RE_EDS/` | s-parameter | ❌ Not implemented |
| **2D EDS-RE** | `replica_exchangers/2D_S_Eoff_RE_EDS/` | s+E_off | ❌ Not implemented |

**Architecture**: Master-slave MPI pattern with replica graph control

**Status**: ❌ **No replica exchange in gromos-rs**

---

## 9. Virtual Atoms & Coarse-Graining

### 9.1 Virtual Sites
| Feature | File | Purpose | Status |
|---------|------|---------|--------|
| **Virtual Atoms** | `virtualatoms/` | Dummy atoms | ❌ Not implemented |
| **Force Propagation** | `virtualatoms/propagate_forces.cc` | Back-distribute forces | ❌ Not implemented |

### 9.2 Coarse-Graining
- Martini-style coarse-grained models
- 4-to-1 mapping (4 atoms → 1 bead)
- Special CG bond/angle potentials

**Status**: ⚠️ **Partial** (basic CG in tests, no full implementation)

---

## 10. Analysis & Properties

### 10.1 Built-in Calculations
| Property | Location | Purpose | Status |
|----------|----------|---------|--------|
| **Energies** | All modules | E_bonded, E_nonbonded, etc. | ⚠️ Partial |
| **Virial/Pressure** | `pressure/pressure_calculation.cc` | Pressure tensor | ❌ Not implemented |
| **Temperature** | `temperature/temperature_calculation.cc` | Kinetic energy | ❌ Not implemented |
| **Averages** | `configuration/average.h` | Running averages | ❌ Not implemented |

### 10.2 gromos++ Tools (separate)
The `gromosPlusPlus` repository contains 100+ analysis programs:
- `rdf` - Radial distribution function
- `hbond` - Hydrogen bond analysis
- `rmsd` - Structure alignment
- `cluster` - Conformational clustering
- `tser` - Time series analysis
- `ener` - Energy analysis
- etc.

**Status**: ❌ **gromos-rs has no analysis tools** (separate project)

---

## 11. GPU Acceleration

### 11.1 CUDA Support
| Feature | File | Purpose | Status |
|---------|------|---------|--------|
| **CUDA Nonbonded** | `nonbonded/interaction/cuda_nonbonded_set.h` | GPU forces | ❌ Not in gromos-rs |
| **GPU SHAKE** | `constraints/gpu_shake.cc` | GPU constraints | ❌ Not in gromos-rs |
| **GPU SETTLE** | `constraints/gpu_settle.cc` | GPU water | ❌ Not in gromos-rs |
| **CUDA Kernels** | `cukernel/` | Low-level GPU code | ❌ Not in gromos-rs |

**Rust Alternative**: Could use `cudarc` or `wgpu` for GPU acceleration

---

## 12. I/O & File Formats

### 12.1 Input Files
| Format | Extension | Purpose | Rust Status |
|--------|-----------|---------|-------------|
| **Topology** | `.top/.topo` | System structure | ✅ Implemented |
| **Coordinates** | `.cnf/.conf` | Positions/velocities | ✅ Implemented |
| **Parameters** | `.imd` | Simulation parameters | ❌ Not implemented |
| **Perturbation** | `.ptp` | Free energy topology | ❌ Not implemented |
| **Restraints** | `.dat/.spec` | Experimental data | ❌ Not implemented |

### 12.2 Output Files
| Format | Extension | Purpose | Rust Status |
|--------|-----------|---------|-------------|
| **Trajectory** | `.trc/.trj` | Coordinates over time | ❌ Not implemented |
| **Energy** | `.tre` | Energy time series | ❌ Not implemented |
| **Forces** | `.trf` | Force output | ❌ Not implemented |
| **Free Energy** | `.dlg` | dH/dλ values | ❌ Not implemented |
| **Restraints** | `.rsr/.dat` | Restraint violations | ❌ Not implemented |

**Status**: gromos-rs can **read** basic files but cannot **write** any GROMOS outputs

---

## 13. Priority Ranking for Implementation

### Tier 1: Core MD (Essential)
1. ✅ **Leap-frog integrator** - Done
2. ✅ **Nonbonded LJ** - Done
3. ❌ **Bonded forces** - CRITICAL MISSING
4. ❌ **SHAKE constraints** - CRITICAL MISSING
5. ❌ **SETTLE (water)** - CRITICAL MISSING
6. ❌ **Thermostat** - CRITICAL MISSING
7. ❌ **Barostat** - CRITICAL MISSING
8. ❌ **Long-range electrostatics** - CRITICAL MISSING

**Estimate**: 4-8 weeks for Tier 1 completion

### Tier 2: Enhanced Methods (Important)
1. ❌ **Steepest descent minimization**
2. ❌ **Stochastic dynamics (Langevin)**
3. ❌ **Distance restraints**
4. ❌ **Position restraints**
5. ❌ **Free energy perturbation**
6. ❌ **Replica exchange**
7. ❌ **Trajectory I/O**

**Estimate**: 8-12 weeks for Tier 2

### Tier 3: Advanced Features (Nice to Have)
1. ❌ **QM/MM**
2. ❌ **X-ray/NMR restraints**
3. ❌ **EDS/GaMD/Metadynamics**
4. ❌ **Virtual atoms**
5. ❌ **GPU acceleration**
6. ❌ **Multi-timestep integrator**

**Estimate**: 16-24 weeks for Tier 3

### Tier 4: Specialized (Optional)
1. ❌ **Monte Carlo**
2. ❌ **NEMD**
3. ❌ **Order parameter restraints**
4. ❌ **Symmetry restraints**
5. ❌ **Perturbed soft-core interactions**

**Estimate**: 12+ weeks for Tier 4

---

## 14. Rust Implementation Strategy

### Phase 1: Complete Core MD (Tier 1)
**Timeline**: 2-3 months
**Deliverable**: Functional MD engine matching basic GROMOS capabilities

**Key Tasks**:
```rust
// Bonded forces
pub mod bonded {
    pub fn calculate_bond_forces(...) -> ForceEnergy;
    pub fn calculate_angle_forces(...) -> ForceEnergy;
    pub fn calculate_dihedral_forces(...) -> ForceEnergy;
}

// Constraints
pub mod constraints {
    pub struct Shake { ... }
    pub struct Settle { ... }
    impl Shake {
        pub fn apply(&mut self, ...) -> Result<()>;
    }
}

// Thermostats
pub mod thermostat {
    pub enum ThermostatType {
        Berendsen { tau: f64 },
        NoseHoover { q: f64 },
    }
    pub trait Thermostat {
        fn apply(&mut self, conf: &mut Configuration, ...);
    }
}

// Long-range electrostatics
pub mod ewald {
    pub fn calculate_ewald(...) -> ForceEnergy;
    pub fn calculate_pme(...) -> ForceEnergy;
}
```

### Phase 2: Enhanced Sampling (Tier 2)
**Timeline**: 2-3 months
**Deliverable**: Free energy, restraints, replica exchange

**Key Features**:
- λ-perturbation for free energy
- Distance/angle/dihedral restraints
- Basic replica exchange (MPI)
- Trajectory output

### Phase 3: Advanced Methods (Tier 3)
**Timeline**: 4-6 months
**Deliverable**: QM/MM, GPU, advanced sampling

**Key Features**:
- QM/MM interface (call external QM)
- GPU acceleration (cudarc)
- EDS, GaMD, metadynamics
- X-ray/NMR refinement

---

## 15. Code Comparison: C++ vs Rust

### Example: Distance Restraint

**C++ GROMOS**:
```cpp
// distance_restraint_interaction.cc
for(auto& restraint : restraints) {
    double r = (conf.pos(restraint.i) - conf.pos(restraint.j)).abs();
    double delta = r - restraint.r0;
    double energy = 0.5 * restraint.k * delta * delta;

    math::Vec f = restraint.k * delta *
                  (conf.pos(restraint.j) - conf.pos(restraint.i)) / r;
    conf.force(restraint.i) += f;
    conf.force(restraint.j) -= f;
}
```

**Rust Implementation**:
```rust
// gromos-rs approach
pub struct DistanceRestraint {
    i: usize,
    j: usize,
    r0: f64,
    k: f64,
}

impl DistanceRestraint {
    pub fn calculate(&self, conf: &mut Configuration) -> f64 {
        let r_vec = conf.pos[self.j] - conf.pos[self.i];
        let r = r_vec.length() as f64;
        let delta = r - self.r0;
        let energy = 0.5 * self.k * delta * delta;

        let f_scalar = self.k * delta / r;
        let force = r_vec * (f_scalar as f32);

        conf.force[self.i] += force;
        conf.force[self.j] -= force;

        energy
    }
}

// Parallel version with Rayon
impl RestraintList {
    pub fn calculate_parallel(&self, conf: &mut Configuration) -> f64 {
        use rayon::prelude::*;

        // Thread-safe force accumulation
        let (total_energy, forces) = self.restraints
            .par_iter()
            .map(|r| r.calculate_force_only(conf))
            .reduce(
                || (0.0, vec![Vec3::ZERO; conf.pos.len()]),
                |(e1, f1), (e2, f2)| (e1 + e2, combine_forces(f1, f2))
            );

        conf.apply_forces(&forces);
        total_energy
    }
}
```

**Rust Advantages**:
- Memory safety (no segfaults)
- Better parallelism (Rayon)
- Clearer ownership semantics
- More composable design

---

## Summary

**GROMOS is MASSIVE**:
- 13 integration algorithms
- 23+ special interactions
- 9 constraint algorithms
- QM/MM hybrid methods
- Replica exchange
- Structure refinement
- Free energy calculations

**gromos-rs Status**:
- ✅ **10% complete** (basic MD only)
- ❌ **90% missing** (most advanced features)

**Recommendation**:
1. **Focus on Tier 1** (core MD) first - 2-3 months
2. Add Tier 2 (enhanced sampling) - another 2-3 months
3. Tier 3+ are "nice to have" for specialized applications

**Timeline to Feature Parity**: ~12-18 months of focused development

**But**: Even with just Tier 1 complete, gromos-rs will be useful for basic MD simulations!
