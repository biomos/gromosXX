# GROMOS Features Catalog - TODO Tracker

**Last Updated**: 2025-11-12
**Purpose**: Track implementation status for gromos-rs

## Legend
- ✅ **Implemented** - Feature complete and tested
- ⚠️ **Partial** - Basic implementation exists, needs work
- ❌ **Not Done** - Not implemented yet
- 🔨 **Effort** - Estimated implementation time

---

## Overall Progress

**Tier 1 (Core MD)**: ✅ 93.75% (7.5/8 complete)
**Tier 2 (Enhanced)**: ✅ 85.7% (6/7 complete)
**Tier 3 (Advanced)**: ❌ 0% (0/6 complete)
**Total Features**: 34% core functionality implemented

---

## 1. Integration Algorithms

| Algorithm | Status | Effort | Location |
|-----------|--------|--------|----------|
| **Leap-Frog** | ✅ Implemented | - | src/integrator.rs:36 |
| **Velocity Verlet** | ✅ Implemented | - | src/integrator.rs:113 |
| **Stochastic Dynamics** | ✅ Implemented | - | src/integrator.rs:387 |
| **Steepest Descent** | ✅ Implemented | - | src/integrator.rs:185 |
| **Slow Growth (FEP)** | ✅ Implemented | - | src/fep.rs:113 |
| **Scaled Leap-Frog** | ❌ Not Done | 🔨 1-2 days | Force scaling + existing leap-frog |
| **Lattice Shift** | ❌ Not Done | 🔨 1-2 weeks | PBC crossing tracker |
| **Conjugate Gradient** | ❌ Not Done | 🔨 2-4 weeks | Needs line search |
| **Monte Carlo** | ❌ Not Done | 🔨 4-6 weeks | Metropolis sampling |
| **Multigradient** | ❌ Not Done | 🔨 4-6 weeks | Multi-potential interpolation |
| **EDS** | ❌ Not Done | 🔨 6-8 weeks | Multi-state sampling |
| **GaMD** | ❌ Not Done | 🔨 6-8 weeks | Boost potential |
| **Analyze** | ❌ Skip | - | Post-processing tool |

**Progress**: 5/13 implemented (38.5%)

---

## 2. Constraint Algorithms

| Algorithm | Status | Effort | Location |
|-----------|--------|--------|----------|
| **SHAKE** | ✅ Implemented | - | src/algorithm/constraints.rs |
| **M-SHAKE** | ✅ Implemented | - | src/algorithm/constraints.rs |
| **SETTLE** | ✅ Implemented | - | src/algorithm/constraints.rs |
| **LINCS** | ✅ Implemented | - | src/algorithm/constraints.rs |
| **Perturbed SHAKE** | ❌ Not Done | 🔨 2-3 weeks | λ-dependent constraints |
| **Flexible Constraints** | ❌ Not Done | 🔨 2 weeks | Time-dependent |
| **COM Motion Removal** | ❌ Not Done | 🔨 1 week | Remove drift |
| **Angle Constraints** | ❌ Not Done | 🔨 2 weeks | Fix angles |
| **Dihedral Constraints** | ❌ Not Done | 🔨 2 weeks | Fix dihedrals |
| **GPU variants** | ❌ Skip | - | Not needed |

**Progress**: 4/9 implemented (44.4%)

---

## 3. Bonded Interactions

| Interaction | Status | Effort | Location |
|-------------|--------|--------|----------|
| **Quartic Bonds** | ✅ Implemented | - | src/interaction/bonded.rs |
| **Harmonic Bonds** | ✅ Implemented | - | src/interaction/bonded.rs |
| **Angles (cosine)** | ✅ Implemented | - | src/interaction/bonded.rs |
| **Proper Dihedrals** | ✅ Implemented | - | src/interaction/bonded.rs |
| **Improper Dihedrals** | ✅ Implemented | - | src/interaction/bonded.rs |
| **Perturbed Terms (FEP)** | ⚠️ Partial | 🔨 2-3 weeks | src/fep.rs (framework done) |
| **Soft-core FEP** | ⚠️ Partial | 🔨 2-3 weeks | src/fep.rs (needs integration) |
| **Harmonic Angles** | ✅ Implemented | - | src/interaction/bonded.rs |
| **CG Bonds** | ✅ Implemented | - | src/interaction/bonded.rs |
| **New Dihedrals** | ✅ Implemented | - | src/interaction/bonded.rs |
| **Cross-Dihedrals** | ✅ Implemented | - | src/interaction/bonded.rs |

**Progress**: 9/11 core terms (81.8%), FEP framework 80% done

---

## 4. Nonbonded Interactions

| Feature | Status | Effort | Location |
|---------|--------|--------|----------|
| **Lennard-Jones** | ✅ Implemented | - | src/interaction/nonbonded.rs |
| **Coulomb (cutoff)** | ✅ Implemented | - | src/interaction/nonbonded.rs |
| **Reaction Field** | ✅ Implemented | - | src/interaction/electrostatics.rs |
| **Pairlist (standard)** | ✅ Implemented | - | src/pairlist.rs |
| **Grid Cell Pairlist** | ✅ Implemented | - | src/pairlist.rs (O(N) spatial decomposition) |
| **Ewald Summation** | ✅ Implemented | - | src/interaction/electrostatics.rs (via PME) |
| **PME** | ✅ Implemented | - | src/interaction/electrostatics.rs (FFT-based Ewald) |
| **P3M** | ❌ Not Done | 🔨 6-8 weeks | Particle-mesh method |

**Progress**: 7/8 (87.5%) - Full long-range electrostatics implemented!
**Note**: RF is GROMOS' traditional long-range method (suitable for most applications)

---

## 5. Special Interactions & Restraints

| Feature | Status | Effort | Category |
|---------|--------|--------|----------|
| **Distance Restraints** | ✅ Implemented | - | src/interaction/restraints.rs |
| **Position Restraints** | ✅ Implemented | - | src/interaction/restraints.rs |
| **Angle Restraints** | ❌ Not Done | 🔨 1-2 weeks | NMR refinement |
| **Dihedral Restraints** | ❌ Not Done | 🔨 1-2 weeks | NMR refinement |
| **J-value Restraints** | ❌ Not Done | 🔨 2-3 weeks | NMR coupling |
| **RDC Restraints** | ❌ Not Done | 🔨 3-4 weeks | NMR orientation |
| **Order Parameter** | ❌ Not Done | 🔨 2-3 weeks | Dynamics |
| **X-ray Restraints** | ❌ Not Done | 🔨 6-8 weeks | Crystallography |
| **Local Elevation** | ❌ Not Done | 🔨 4-6 weeks | Metadynamics |
| **Distance Field** | ❌ Not Done | 🔨 3-4 weeks | Biasing potential |
| **Electric Field** | ❌ Not Done | 🔨 1-2 weeks | External E-field |
| **NEMD** | ❌ Not Done | 🔨 3-4 weeks | Non-equilibrium |
| **Symmetry Restraints** | ❌ Not Done | 🔨 2-3 weeks | Enforce symmetry |

**Progress**: 2/13 implemented (15.4%)

---

## 6. Thermostats & Barostats

| Feature | Status | Location |
|---------|--------|----------|
| **Berendsen Thermostat** | ✅ Implemented | src/algorithm/thermostats.rs |
| **Nosé-Hoover Thermostat** | ✅ Implemented | src/algorithm/thermostats.rs |
| **Andersen Thermostat** | ✅ Implemented | src/algorithm/thermostats.rs |
| **Berendsen Barostat** | ✅ Implemented | src/algorithm/barostats.rs |
| **Parrinello-Rahman Barostat** | ✅ Implemented | src/algorithm/barostats.rs |

**Progress**: 5/5 implemented (100%)

---

## 7. Advanced Features (Not Yet Implemented)

| Feature Category | Status | Effort | Priority |
|------------------|--------|--------|----------|
| **QM/MM** | ❌ Not Done | 🔨 12+ weeks | Tier 3 |
| **Replica Exchange (REMD)** | ❌ Not Done | 🔨 4-6 weeks | Tier 2 |
| **Virtual Atoms** | ❌ Not Done | 🔨 3-4 weeks | Tier 3 |
| **Coarse-Graining** | ⚠️ Partial | 🔨 4-6 weeks | Tier 3 |
| **GPU Acceleration** | ❌ Not Done | 🔨 8-12 weeks | Tier 3 |
| **Analysis Tools** | ❌ Not Done | Variable | Separate project |

**Note**: These are advanced/specialized features for future implementation

---

## 8. I/O & File Formats

| Format | Extension | Status | Location |
|--------|-----------|--------|----------|
| **Topology (read)** | `.top` | ✅ Implemented | src/io/topology.rs |
| **Coordinates (read)** | `.cnf` | ✅ Implemented | src/io/coordinate.rs |
| **Trajectory (write)** | `.trc` | ✅ Implemented | src/io/trajectory.rs |
| **Energy (write)** | `.tre` | ✅ Implemented | src/io/energy.rs |
| **Forces (write)** | `.trf` | ✅ Implemented | src/io/force.rs |
| **Parameters** | `.imd` | ✅ Implemented | src/io/imd.rs (426 lines, full parser) |
| **Perturbation** | `.ptp` | ❌ Not Done | FEP topology |
| **Free Energy** | `.dlg` | ❌ Not Done | dH/dλ output |
| **Restraints** | `.dat/.spec` | ❌ Not Done | Experimental data |

**Progress**: 6/9 formats (66.7%) - Core I/O complete, FEP files remaining

---

## 9. Implementation Priorities

### ✅ Tier 1: Core MD (93.75% Complete)
**Status**: Nearly complete, production-ready for basic MD

| Feature | Status |
|---------|--------|
| Integrators (Leap-frog, Verlet) | ✅ Done |
| Bonded forces (all standard terms) | ✅ Done |
| Nonbonded (LJ + Coulomb) | ✅ Done |
| Constraints (SHAKE, SETTLE) | ✅ Done |
| Thermostats (3 types) | ✅ Done |
| Barostats (2 types) | ✅ Done |
| Reaction Field electrostatics | ⚠️ Complete this (1-2 weeks) |

### ⚠️ Tier 2: Enhanced Methods (85.7% Complete)
**Status**: Most features done, 1 major item remaining

| Feature | Status |
|---------|--------|
| Energy minimization | ✅ Done |
| Stochastic dynamics | ✅ Done |
| Distance/Position restraints | ✅ Done |
| Free energy perturbation | ✅ Done (framework) |
| Trajectory I/O | ✅ Done |
| **Replica Exchange** | ❌ TODO (4-6 weeks) |

### ❌ Tier 3: Advanced Features (0% Complete)
**Next Steps**: Pick based on scientific needs

| Feature | Effort | Use Case |
|---------|--------|----------|
| PME long-range | 6-8 weeks | Charged systems |
| Conjugate Gradient | 2-4 weeks | Fast minimization |
| QM/MM | 12+ weeks | Reactive systems |
| EDS/GaMD | 6-8 weeks each | Enhanced sampling |
| Virtual atoms | 3-4 weeks | Special topologies |

### Quick Wins (1-2 weeks each)
1. **Scaled Leap-Frog** (1-2 days) - Multiple time-stepping
2. **Lattice Shift** (1-2 weeks) - FEP with long-range
3. **Grid Cell Pairlist** (2-3 weeks) - Performance boost

---

## 10. Current Status Summary

### What Works Now (Production-Ready)
- ✅ **Core MD Engine**: NVE, NVT, NPT ensembles
- ✅ **Force Field**: All standard bonded/nonbonded terms
- ✅ **Constraints**: SHAKE, M-SHAKE, SETTLE
- ✅ **Temperature/Pressure Control**: 3 thermostats, 2 barostats
- ✅ **Minimization**: Steepest descent
- ✅ **Stochastic Dynamics**: Langevin integrator
- ✅ **Restraints**: Position and distance restraints
- ✅ **Free Energy**: FEP framework with soft-core potentials
- ✅ **I/O**: Read TOP/CNF, write TRC/TRE/TRF trajectories

### What's Missing (Top Priorities)
1. **Complete Reaction Field** (1-2 weeks) - Finish Tier 1
2. **Replica Exchange** (4-6 weeks) - Complete Tier 2
3. **Quick wins**: Scaled Leap-Frog (1-2 days), Grid Pairlist (2-3 weeks)

### Advanced Features (Future Work)
- PME long-range electrostatics (6-8 weeks)
- Conjugate Gradient minimization (2-4 weeks)
- QM/MM, EDS, GaMD (12+ weeks each)
- Virtual atoms, GPU acceleration

### Tutorial Tools & Binaries
- ✅ **Implemented**: pdb2g96, com_top, check_top, mk_script, md, ene_ana, rmsd (7+ tools)
  - **md** (946 lines): Full MD simulation engine with CLI
  - **ene_ana**: Energy trajectory analysis with statistics
  - **rmsd**: RMSD calculation with optional rotational fitting
- ❌ **Needed**: sim_box (box generation utility)

---

## Next Actions

### To Complete Tier 1 (1-2 weeks)
- Finish Reaction Field implementation and testing

### To Complete Tier 2 (4-6 weeks)
- Implement Replica Exchange (T-REMD, H-REMD)

### Quick Performance Wins (1-3 weeks each)
- Scaled Leap-Frog integrator
- Grid cell pairlist algorithm
- Lattice shift tracking
