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

**Tier 1 (Core MD)**: ✅ **100%** (8/8 complete) - PRODUCTION READY!
**Tier 2 (Enhanced)**: ✅ **85.7%** (6/7 complete)
**Tier 3 (Advanced)**: ⚠️ **16.7%** (1/6 complete) - PME implemented
**Total Features**: **~60%** of core functionality implemented

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
| **Perturbed Bonds (FEP)** | ✅ Implemented | - | src/interaction/bonded.rs (with λ derivatives) |
| **Perturbed Angles (FEP)** | ❌ Not Done | 🔨 1-2 weeks | Similar to perturbed bonds |
| **Perturbed Dihedrals (FEP)** | ❌ Not Done | 🔨 1-2 weeks | Similar to perturbed bonds |
| **Soft-core FEP** | ⚠️ Partial | 🔨 1-2 weeks | src/fep.rs (framework ready, needs nonbonded integration) |
| **Harmonic Angles** | ✅ Implemented | - | src/interaction/bonded.rs |
| **CG Bonds** | ✅ Implemented | - | src/interaction/bonded.rs |
| **New Dihedrals** | ✅ Implemented | - | src/interaction/bonded.rs |
| **Cross-Dihedrals** | ✅ Implemented | - | src/interaction/bonded.rs |

**Progress**: 9/11 core terms (81.8%), FEP: perturbed bonds ✅, framework 85% done

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

| Format | Extension | Status | Location | Notes |
|--------|-----------|--------|----------|-------|
| **Topology (read)** | `.top` | ✅ Implemented | src/io/topology.rs | Input for simulations |
| **Coordinates (read)** | `.cnf` | ✅ Implemented | src/io/coordinate.rs | Initial positions |
| **Parameters (read)** | `.imd` | ✅ Implemented | src/io/imd.rs | 426 lines, full parser |
| **Trajectory (write)** | `.trc` | ✅ Implemented | src/io/trajectory.rs | Positions/velocities/forces |
| **Energy (write)** | `.tre` | ✅ Implemented | src/io/energy.rs | Energy components |
| **Forces (write)** | `.trf` | ✅ Implemented | src/io/force.rs | Detailed force breakdown |
| **Perturbation (write)** | `.ptp` | ✅ Implemented | src/io/ptp.rs | FEP topology writer (make_pt_top tool) |
| **Free Energy** | `.dlg` | ❌ Not Done | dH/dλ output writer | TI analysis |
| **Restraints** | `.dat/.spec` | ❌ Not Done | NMR/experimental data | Future |

**Progress**: 7/9 formats (77.8%) - FEP I/O nearly complete!

**Note**: Trajectory *reading* is handled by GROMOS++ (111 analysis tools). gromos-rs focuses on *writing* trajectories efficiently.

---

## 9. Implementation Priorities

### ✅ Tier 1: Core MD (100% Complete) - PRODUCTION READY!
**Status**: ✅ **COMPLETE** - Fully functional MD engine ready for simulations

| Feature | Status |
|---------|--------|
| Integrators (Leap-frog, Verlet, SD) | ✅ Done |
| Bonded forces (all 11 standard terms) | ✅ Done |
| Nonbonded (LJ + Coulomb + RF) | ✅ Done |
| Constraints (SHAKE, M-SHAKE, SETTLE, LINCS) | ✅ Done |
| Thermostats (Berendsen, Nosé-Hoover, Andersen) | ✅ Done |
| Barostats (Berendsen, Parrinello-Rahman) | ✅ Done |
| Reaction Field electrostatics | ✅ Done |
| Grid Cell Pairlist (O(N) performance) | ✅ Done |

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

### ⚠️ Tier 3: Advanced Features (16.7% Complete)
**Next Steps**: Pick based on scientific needs

| Feature | Status | Effort | Use Case |
|---------|--------|--------|----------|
| **PME long-range** | ✅ Done | - | Charged systems (fully functional!) |
| **Ewald Summation** | ✅ Done | - | Periodic electrostatics |
| Conjugate Gradient | ❌ TODO | 🔨 2-4 weeks | Fast minimization |
| QM/MM | ❌ TODO | 🔨 12+ weeks | Reactive systems |
| EDS/GaMD | ❌ TODO | 🔨 6-8 weeks each | Enhanced sampling |
| Virtual atoms | ❌ TODO | 🔨 3-4 weeks | Special topologies |

### Quick Wins Remaining (1-2 weeks each)
1. **Scaled Leap-Frog** (1-2 days) - Multiple time-stepping
2. **Lattice Shift** (1-2 weeks) - FEP with long-range
3. **Angle/Dihedral Restraints** (1-2 weeks each) - NMR refinement

---

## 10. Current Status Summary

### What Works Now (Production-Ready) ✅

**Core MD Engine**:
- ✅ All ensembles: NVE, NVT, NPT
- ✅ Integrators: Leap-frog, Velocity Verlet, Stochastic Dynamics
- ✅ Minimization: Steepest Descent

**Force Field (Complete)**:
- ✅ Bonded: 11/11 terms (quartic/harmonic bonds, cosine/harmonic angles, proper/improper dihedrals, cross-dihedrals, CG bonds, new dihedrals)
- ✅ Nonbonded: Lennard-Jones, Coulomb (cutoff)
- ✅ Long-range: Reaction Field, PME/Ewald Summation
- ✅ Pairlist: Standard + Grid Cell (O(N) performance)

**Constraints**:
- ✅ SHAKE, M-SHAKE, SETTLE, LINCS

**Thermostats & Barostats**:
- ✅ Thermostats: Berendsen, Nosé-Hoover, Andersen
- ✅ Barostats: Berendsen, Parrinello-Rahman

**Free Energy Perturbation**:
- ✅ FEP framework with lambda control
- ✅ Perturbed bonds with lambda derivatives
- ✅ Soft-core potentials (framework ready)

**Restraints**:
- ✅ Position and distance restraints

**I/O (Complete)**:
- ✅ Read: .top (topology), .cnf (coordinates), .imd (parameters)
- ✅ Write: .trc (trajectory), .tre (energy), .trf (forces)

**Analysis Tools**:
- ✅ md (full MD simulation binary, 946 lines)
- ✅ ene_ana (energy analysis)
- ✅ rmsd (RMSD calculator)
- ✅ pdb2g96, com_top, check_top, mk_script, and more

### What's Missing (Top Priorities)
1. **Replica Exchange (T-REMD, H-REMD)** (4-6 weeks) - Complete Tier 2
2. **FEP I/O (.ptp, .dlg readers/writers)** (2-3 weeks) - FEP simulations
3. **Quick wins**: Scaled Leap-Frog (1-2 days), Lattice Shift (1-2 weeks)

### Advanced Features (Future Work)
- Conjugate Gradient minimization (2-4 weeks)
- QM/MM, EDS, GaMD (12+ weeks each)
- Virtual atoms, GPU acceleration
- P3M (alternative to PME, 6-8 weeks)
- NMR restraints (J-value, RDC, 2-4 weeks each)

### Tutorial Tools & Binaries

**gromos-rs Binaries** (Simulation & Conversion):
- ✅ **md** (946 lines): Full MD simulation engine with CLI
- ✅ **make_pt_top** (265 lines): Generate .ptp perturbation topologies for FEP
- ✅ **pdb2g96**: Convert PDB to GROMOS format
- ✅ **com_top**: Combine topology files
- ✅ **check_top**: Validate topology
- ✅ **mk_script**: Generate simulation scripts
- ✅ **frameout**, **trs_ana**, **diffus**, **hbond**, **rdf**, **rmsf**: Basic analysis tools
- ❌ **sim_box**: Box generation utility (solvation, nice to have)

**Analysis Strategy**:
Use **GROMOS++** for advanced analysis (111 battle-tested tools):
- `ener_ana`, `rmsd`, `hbond`, `cluster`, `rdf`, `sasa`, `bar`, etc.
- gromos-rs writes .trc/.tre/.trf → GROMOS++ reads and analyzes
- See `GROMOS_PLUSPLUS_INTEGRATION.md` for details
- **Don't reimplement** - leverage existing tools!

---

## Next Actions

### ✅ Tier 1 Complete!
All Tier 1 features are now implemented and production-ready!

### To Complete Tier 2 (4-6 weeks)
- Implement Replica Exchange (T-REMD, H-REMD)

### Quick Wins (1-3 weeks each)
- ✅ ~~Grid cell pairlist algorithm~~ - **DONE**
- ✅ ~~PME/Ewald Summation~~ - **DONE**
- Scaled Leap-Frog integrator (1-2 days)
- Lattice shift tracking (1-2 weeks)
- NMR restraints (angle/dihedral, 1-2 weeks each)

### FEP Completion (2-3 weeks)
- .ptp perturbation topology reader
- .dlg free energy output writer
- Complete perturbed angles and dihedrals (bonds already done!)

---

## 11. GROMOS++ Integration Strategy

**Philosophy**: Don't reimplement - integrate!

### Division of Labor
| Component | Responsibility | Status |
|-----------|---------------|--------|
| **gromos-rs** | Simulation engine (write trajectories) | ✅ Production ready |
| **GROMOS++** | Analysis tools (read trajectories) | ✅ Use existing 111 tools |

### Workflow
```
gromos-rs (Rust)        GROMOS++ (C++)
    md binary     -->    .trc/.tre/.trf    -->    ener_ana, rmsd, hbond
   (simulate)           (trajectories)             (analyze)
```

### What gromos-rs Does
- ✅ Run MD simulations (NVE/NVT/NPT)
- ✅ Write trajectory files (.trc, .tre, .trf)
- ✅ Topology conversion (pdb2g96, com_top)
- ✅ Basic utilities (check_top, mk_script)

### What GROMOS++ Does (Use Existing Tools)
- Energy analysis: `ener_ana`, `int_ener`, `dg_ener`
- Structural: `rmsd`, `rmsf`, `rgyr`, `dssp`, `sasa`
- Interactions: `hbond`, `rdf`, `ion`, `close_pair`
- Clustering: `cluster`, `follow`
- Free energy: `bar`, `ext_ti_ana`
- **111 tools total** - battle-tested, maintained

### Building GROMOS++
```bash
cd gromosPlusPlus/gromos++
./Config.sh && ./configure && make -j$(nproc)
# Binaries in: programs/
```

### Why This Strategy?
1. ✅ Avoid duplication (20+ years of development)
2. ✅ Proven, tested code
3. ✅ Focus gromos-rs on simulation
4. ✅ User familiarity (existing workflows)
5. ✅ Less maintenance burden

See `GROMOS_PLUSPLUS_INTEGRATION.md` for complete details.
