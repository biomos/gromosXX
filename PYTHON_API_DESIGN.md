# GROMOS-RS Codebase Analysis & Python API Recommendations

## Executive Summary

The **gromos-rs** project is a high-performance Rust implementation of GROMOS (GROningen MOlecular Simulation), a molecular dynamics engine for biomolecular simulations. With ~7,500 lines of Rust code, it provides:

- **Core MD Engine**: Fully functional simulation with multiple integrators
- **Advanced Sampling**: GAMD, EDS, and REMD capabilities
- **Advanced Features**: FEP (Free Energy Perturbation), constraints, restraints
- **Analysis Tools**: 23 command-line analysis tools (RMSD, RMSF, RMSF, RDF, energy analysis, etc.)
- **Performance**: SIMD acceleration, parallelization with Rayon, grid-cell pairlist (O(N))

**Crate Status**: Production-ready core MD (~60% of total feature catalog)

---

## Architecture Overview

```
gromos-rs
├── Core Systems
│   ├── math              - SIMD vectors (Vec3, Mat3) with glam+wide
│   ├── topology          - Molecular structure (atoms, bonds, angles, etc.)
│   ├── configuration     - System state (positions, velocities, forces, energies)
│   ├── pairlist          - Neighbor lists (standard + grid-cell O(N))
│   └── interaction       - Force calculations
│       ├── bonded        - Bonds, angles, dihedrals
│       ├── nonbonded     - LJ + Coulomb with soft-core FEP
│       ├── electrostatics - CRF, PME, Ewald summation
│       └── restraints    - Distance, position, angle, dihedral restraints
│
├── Simulation Control
│   ├── integrator        - Leap-Frog, Velocity Verlet, Stochastic Dynamics, Steepest Descent
│   ├── algorithm
│   │   ├── thermostats   - Berendsen, Nosé-Hoover, Andersen
│   │   ├── barostats     - Berendsen, Parrinello-Rahman
│   │   └── constraints   - SHAKE, M-SHAKE, SETTLE, LINCS
│   └── fep               - Free Energy Perturbation (Thermodynamic Integration, Slow Growth)
│
├── Enhanced Sampling
│   ├── gamd              - Gaussian Accelerated MD (3 search modes, dual boost)
│   ├── eds               - Enveloping Distribution Sampling (multi-state, adaptive)
│   ├── remd              - Replica Exchange MD (temperature/lambda/2D exchange)
│   └── replica           - Individual replica management
│
├── I/O System
│   └── io
│       ├── topology      - .top file reading
│       ├── coordinate    - .cnf file reading
│       ├── imd           - .imd parameter file parsing
│       ├── trajectory    - .trc trajectory writing
│       ├── energy        - .tre energy frame writing
│       ├── force         - .trf force output
│       ├── ptp           - .ptp perturbation topology (FEP)
│       ├── dlg           - .dlg lambda derivatives (TI)
│       └── output        - GAMD/EDS statistics writers
│
├── Utilities
│   ├── selection         - Atom selection syntax (all, ranges, residues, names)
│   ├── validation        - Input validation
│   ├── logging           - Structured logging + diagnostics
│   └── (23 analysis tools as separate binaries)
```

**Code Statistics**:
- Total: ~7,500 lines of Rust
- Largest modules: bonded.rs (2,624 lines), electrostatics.rs (1,137 lines), eds.rs (1,036 lines)
- Compiled outputs: static lib + cdylib + rlib (Python binding-ready)

---

## 1. Main Data Structures for Python Exposure

### 1.1 Core Molecular System

#### `Configuration` - System State Container
**File**: `src/configuration.rs`
**Key Members**:
- `current()` / `current_mut()` - Get current state
- `old()` / `old_mut()` - Get previous state (double buffering)
- `exchange_state()` - Zero-cost pointer swap

#### `State` - Single Time Point
**Key Members**:
```rust
pub struct State {
    pub pos: Vec<Vec3>,              // Atomic positions [nm]
    pub vel: Vec<Vec3>,              // Velocities [nm/ps]
    pub force: Vec<Vec3>,            // Forces [kJ/(mol·nm)]
    pub constraint_force: Vec<Vec3>, // Constraint forces
    pub box_config: Box,             // Simulation box
    pub virial_tensor: Mat3,         // For pressure calculation
    pub energies: Energy,            // Energy components
}
```
**Python Usage**: Access/modify/read simulation state

#### `Energy` - Detailed Energy Breakdown
**Key Members**:
```rust
pub struct Energy {
    // Totals
    pub kinetic_total: f64,
    pub potential_total: f64,
    
    // Bonded energies
    pub bond_total: f64,
    pub angle_total: f64,
    pub dihedral_total: f64,
    pub improper_total: f64,
    pub cross_dihedral_total: f64,
    
    // Nonbonded energies
    pub lj_total: f64,               // Lennard-Jones
    pub crf_total: f64,              // Coulomb Reaction Field
    pub ls_total: f64,               // Lattice sum (PME/Ewald)
    pub special_total: f64,
    pub sasa_total: f64,
    
    // Per-group energies
    pub kinetic_energy: Vec<f64>,    // [num_temperature_groups]
    pub lj_energy: Vec<Vec<f64>>,    // [ng x ng] matrix
    pub crf_energy: Vec<Vec<f64>>,   // [ng x ng] matrix
    
    pub virial_total: f64,
}
```
**Python Usage**: Query energy components, temperature calculation

#### `Box` - Periodic Boundary Conditions
**Key Members**:
```rust
pub enum BoxType { Vacuum, Rectangular, Triclinic, TruncatedOctahedral }
pub struct Box {
    pub box_type: BoxType,
    pub vectors: Mat3,      // Box matrix columns
    pub inv_vectors: Mat3,  // Inverse for wrapping
}
```
**Methods**: `volume()`, `dimensions()`, `rectangular()`, `triclinic()`

### 1.2 Molecular Structure

#### `Topology` - Complete Molecular System
**File**: `src/topology.rs` (1,297 lines total)
**Key Members**:
```rust
pub struct Topology {
    pub atoms: Vec<Atom>,                           // [n_atoms]
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub dihedrals: Vec<Dihedral>,
    pub improper_dihedrals: Vec<Dihedral>,
    pub cross_dihedrals: Vec<CrossDihedral>,
    
    // FEP support
    pub perturbed_atoms: Vec<usize>,
    pub perturbed_bonds: Vec<PerturbedBond>,
    pub perturbed_angles: Vec<PerturbedAngle>,
    pub perturbed_dihedrals: Vec<PerturbedDihedral>,
    
    // Force field parameters
    pub bond_parameters: Vec<BondParameters>,
    pub angle_parameters: Vec<AngleParameters>,
    pub dihedral_parameters: Vec<DihedralParameters>,
    pub improper_parameters: Vec<ImproperDihedralParameters>,
    pub lj_parameters: Vec<Vec<LJParameters>>,      // [iac x iac] matrix
    pub charge: Vec<f64>,                           // [n_atoms]
    pub mass: Vec<f64>,                             // [n_atoms]
    pub inverse_mass: Vec<f64>,                     // [n_atoms] (1/mass)
    
    // Metadata
    pub num_atoms: usize,
    pub num_molecules: usize,
}
```
**Methods**: `num_atoms()`, `num_molecules()`, parse from file

#### `Atom` - Individual Atom Properties
```rust
pub struct Atom {
    pub name: String,
    pub residue_nr: usize,
    pub residue_name: String,
    pub iac: usize,           // Integer atom code (type)
    pub mass: f64,
    pub charge: f64,
    pub is_perturbed: bool,
    pub is_polarisable: bool,
    pub is_coarse_grained: bool,
}
```

### 1.3 Math Primitives

#### `Vec3` & `Mat3` - SIMD-Accelerated Vectors/Matrices
**From glam library**, re-exported by gromos-rs
```rust
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

pub struct Mat3 {
    pub x_axis: Vec3,
    pub y_axis: Vec3,
    pub z_axis: Vec3,
}
```
**Methods**: `dot()`, `cross()`, `length()`, matrix operations
**Python Usage**: Convert to numpy arrays

---

## 2. Simulation & Integration Capabilities

### 2.1 Integrators (src/integrator.rs)

| Algorithm | Status | Use Case |
|-----------|--------|----------|
| **LeapFrog** | ✅ | Default, best energy conservation |
| **VelocityVerlet** | ✅ | Alternative integrator |
| **StochasticDynamics** | ✅ | Langevin dynamics, implicit solvent |
| **SteepestDescent** | ✅ | Minimization |

**Trait**: `Integrator` with `step(dt, topo, conf)` method
**Python Usage**: Switch integrators dynamically

### 2.2 Temperature & Pressure Control

#### Thermostats (src/algorithm/thermostats.rs)
- **Berendsen**: Simple, fast, weak ensemble
- **Nosé-Hoover**: Proper NVT ensemble
- **Andersen**: Stochastic, random collision events

#### Barostats (src/algorithm/barostats.rs)
- **Berendsen**: Simple, fast scaling
- **Parrinello-Rahman**: Proper NPT ensemble

**Python Usage**: Apply coupling during dynamics

### 2.3 Constraints (src/algorithm/constraints.rs)

| Type | Status | Atoms | Use |
|------|--------|-------|-----|
| **SHAKE** | ✅ | 2 | Bond length |
| **M-SHAKE** | ✅ | 3+ | Multi-atom |
| **SETTLE** | ✅ | 3 | Water geometry |
| **LINCS** | ✅ | 2-3 | Accurate, stable |

**Python Usage**: Enforce bond/angle constraints

---

## 3. Advanced Sampling Methods

### 3.1 Gaussian Accelerated MD (GaMD)

**Files**: `src/gamd.rs` (864 lines), `src/bin/gamd.rs` (532 lines)

**Features**:
- **Three search modes**: NoSearch (production), CmdSearch (statistics collection), GamdSearch (adaptive boost)
- **Boost forms**: DihedralBoost, TotalBoost, DualBoost
- **Threshold types**: LowerBound, UpperBound
- **Online statistics**: Welford's algorithm for mean/variance

**Key Structures**:
```rust
pub enum SearchMode { NoSearch, CmdSearch, GamdSearch }
pub enum BoostForm { DihedralBoost, TotalBoost, DualBoost }
pub enum ThresholdType { LowerBound, UpperBound }

pub struct GamdParameters {
    pub search_mode: SearchMode,
    pub boost_form: BoostForm,
    pub threshold_type: ThresholdType,
    pub dih_stats: GamdStatistics,  // Dihedral energy statistics
    pub tot_stats: GamdStatistics,  // Total energy statistics
    pub sigma0_dih: f64,            // Max std dev (dihedral)
    pub sigma0_tot: f64,            // Max std dev (total)
    pub k0_dih: f64,                // Force constant factor
    pub k0_tot: f64,
}

pub struct GamdStatistics {
    pub v_max: f64,
    pub v_min: f64,
    pub v_mean: f64,
    pub sigma_v: f64,
    pub steps: usize,
    pub is_ready() -> bool          // Can calculate parameters?
}

pub struct GamdRunner {
    // Manages full GaMD simulation
}
```

**Workflow**:
1. Classical MD (CmdSearch) - collect statistics
2. Accelerated search (GamdSearch) - adaptive parameters
3. Production (NoSearch) - fixed parameters

**Python Usage**: Run GaMD simulations with automatic parameter adaptation

### 3.2 Enveloping Distribution Sampling (EDS)

**Files**: `src/eds.rs` (1,036 lines), `src/bin/eds.rs` (503 lines)

**Features**:
- **Multi-state sampling**: 2+ end-states
- **Smooth Hamiltonian**: Reference potential envelopes all states
- **Adaptive EDS (AEDS)**: Automatic offset updates
- **Flexible s-values**: SingleS, MultiS, PairS forms

**Key Structures**:
```rust
pub enum EDSForm { SingleS, MultiS, PairS }

pub struct EDSState {
    pub id: usize,
    pub energy: f64,                // V_i(R)
    pub offset: f64,                // E_i^R
    pub forces: Vec<Vec3>,          // Forces for this state
    pub virial: [[f64; 3]; 3],
    pub visit_count: usize,
    pub visited: bool,
    pub adjusted_energy(&self) -> f64  // V_i - E_i^R
}

pub struct EDSParameters {
    pub form: EDSForm,
    pub num_states: usize,
    pub s_values: Vec<f64>,         // Smoothing parameters
    pub temperature: f64,
    pub states: Vec<EDSState>,
    pub reference_energy: f64,      // V_R
    pub mixed_energy: f64,
    pub current_state_id: usize,
    pub round_trips: usize,
    pub beta(&self) -> f64          // 1/(kT)
}

pub struct AEDSParameters {
    // Adaptive offsets (auto-update)
    pub update_frequency: usize,
    pub offset_target: f64,         // Target probability for each state
}

pub struct EDSRunner {
    // Manages full EDS simulation
}
```

**Method**: Log-sum-exp for numerical stability, Metropolis round-trip analysis

**Python Usage**: Sample across multiple molecular states (ligand binding, conformations)

### 3.3 Replica Exchange MD (REMD)

**Files**: `src/remd.rs` (1,277 lines), `src/bin/remd.rs` (535 lines)

**Features**:
- **Exchange types**: Temperature only, Lambda (FEP), 2D Temp-Lambda
- **Exchange schemes**: Sequential, OddEven, Random
- **Statistics tracking**: Per-pair acceptance rates
- **Parallel replicas**: Using Rayon

**Key Structures**:
```rust
pub enum ExchangeType { Temperature, Lambda, TemperatureLambda }
pub enum ExchangeScheme { Sequential, OddEven, Random }

pub struct ReplicaInfo {
    pub id: ReplicaId,
    pub temperature: f64,
    pub lambda: f64,
    pub dt: f64,
    pub potential_energy: f64,
    pub partner_id: Option<ReplicaId>,
    pub exchange_probability: f64,
    pub exchange_accepted: bool,
    pub run: usize,
    pub volume: f64,
    pub beta(&self) -> f64          // 1/(kT)
}

pub struct Replica {
    pub info: ReplicaInfo,
    pub configuration: Configuration,
    pub lambda_controller: Option<LambdaController>,
    // ... internal state ...
}

pub struct ExchangeStatistics {
    pub attempts: usize,
    pub accepted: usize,
    pub acceptance_rate: f64,
    pub pair_acceptances: Vec<Vec<usize>>,  // [i][j]
    pub pair_attempts: Vec<Vec<usize>>,     // [i][j]
    pub record_attempt(id1, id2, accepted: bool)
    pub pair_acceptance_rate(id1, id2) -> f64
}

pub struct ReplicaController {
    pub replicas: Vec<Replica>,
    pub exchange_type: ExchangeType,
    pub exchange_scheme: ExchangeScheme,
    pub exchange_interval: usize,
    pub equilibration_steps: usize,
    pub statistics: ExchangeStatistics,
}
```

**Workflow**: 
1. Run each replica in parallel
2. At intervals, propose exchanges between pair of replicas
3. Accept/reject with Metropolis criterion
4. Track acceptance statistics for each pair

**Python Usage**: Enhanced sampling via parallel tempering or FEP

---

## 4. Free Energy Perturbation (FEP)

**File**: `src/fep.rs` (687 lines)

**Features**:
- **Lambda-dependent potentials**: Interpolate between state A and B
- **Soft-core potentials**: Prevent singularities in LJ/Coulomb
- **Individual lambda values**: Per interaction type (bonds, LJ, CRF, etc.)
- **Thermodynamic Integration**: Collect dH/dλ for free energy calculation
- **Slow Growth**: Continuous lambda evolution

**Key Structures**:
```rust
pub struct LambdaController {
    pub lambda: f64,                // Current coupling parameter [0,1]
    pub dlambda_dt: f64,            // Evolution rate (slow growth)
    pub lambda_exponent: i32,       // Power law coupling
    pub interaction_lambdas: InteractionLambdas,
    pub update_lambda(dt: f64)      // Time evolution
}

pub struct InteractionLambdas {
    pub bond: f64,
    pub angle: f64,
    pub dihedral: f64,
    pub improper: f64,
    pub lj: f64,
    pub lj_softness: f64,           // Soft-core strength
    pub crf: f64,
    pub crf_softness: f64,          // Soft-core strength
    pub mass: f64,
}

pub struct SoftCoreParameters {
    pub alpha: f64,                 // Soft-core exponent (1-2)
    pub lambda_power: i32,          // λ^n
    pub sc_sigma: f64,              // Min distance cutoff
}

pub struct FreeEnergyDerivatives {
    pub dH_dlambda: f64,            // ⟨∂H/∂λ⟩
    pub dH_dlambda_dih: f64,        // Dihedral contribution
    pub dH_dlambda_nonbond: f64,    // Nonbonded contribution
}

pub struct PerturbedTopology {
    pub n_perturbed: usize,
    pub lambda_derivative: FreeEnergyDerivatives,
}
```

**Python Usage**: Compute free energies using TI or slow-growth FEP

---

## 5. Interaction Calculations

### 5.1 Bonded Interactions (src/interaction/bonded.rs - 2,624 lines)

**Supported Terms**:
- **Bonds**: Quartic, Harmonic
- **Angles**: Harmonic, Harmonic (with cosine option)
- **Dihedrals**: Proper (periodic), Improper
- **Cross-dihedrals**: 8-atom terms
- **All with FEP support**: State A/B parameters with soft-core

**Functions**:
```rust
pub fn calculate_bond_forces_quartic(...) -> (f64, f64, f64)  // (energy, dV/dr_0, dV/dr_1)
pub fn calculate_bond_forces_harmonic(...)
pub fn calculate_angle_forces(...)
pub fn calculate_dihedral_forces(...)
pub fn calculate_improper_dihedral_forces(...)
pub fn calculate_bonded_forces(topo, conf) -> f64  // Returns total bonded energy
pub fn calculate_perturbed_*(...) -> (f64_A, f64_B, dH_dlambda)  // FEP versions
```

**Python Usage**: Access/debug individual force components

### 5.2 Nonbonded Interactions (src/interaction/nonbonded.rs - 907 lines)

**Methods**:
- **Lennard-Jones**: C6/r^6 - C12/r^12 (standard + soft-core)
- **Coulomb Reaction Field**: CRF long-range (GROMOS standard)
- **Cutoff handling**: Shift force/energy at cutoff

**Key Functions**:
```rust
pub fn lj_crf_innerloop(...)    // Single-threaded
pub fn lj_crf_innerloop_parallel(...)  // Rayon parallelized
pub struct ForceStorage {
    pub force: Vec<Vec3>,
    pub virial: [[f64; 3]; 3],
}
```

**Features**:
- Grid-cell pairlist (O(N) scaling)
- Parallel force calculation
- Energy group-based accounting
- FEP with soft-core

### 5.3 Electrostatics (src/interaction/electrostatics.rs - 1,137 lines)

**Methods**:
- **Reaction Field (CRF)**: GROMOS standard, fast
- **Particle Mesh Ewald (PME)**: Long-range, FFT-based
- **Ewald Summation**: Full implementation

**Key Parameters**:
```rust
pub struct ReactionFieldParameters {
    pub cutoff: f64,
    pub epsilon_inside: f64,         // Dielectric inside cutoff
    pub epsilon_outside: f64,        // Dielectric outside cutoff (infinity for CRF)
    pub kappa: f64,                  // Screening parameter
}

pub struct PMEParameters {
    pub cutoff: f64,
    pub pme_order: usize,            // B-spline order
    pub pme_grid_x: usize,
    pub pme_grid_y: usize,
    pub pme_grid_z: usize,
    pub ew_coeff: f64,               // Ewald coefficient
}
```

### 5.4 Restraints (src/interaction/restraints.rs - 858 lines)

**Types**:
```rust
pub struct DistanceRestraint {
    pub i: usize,
    pub j: usize,
    pub r0: f64,                     // Reference distance
    pub k: f64,                      // Force constant
}

pub struct PositionRestraint {
    pub atom_idx: usize,
    pub reference_pos: Vec3,
    pub k: f64,
}

pub struct AngleRestraint {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub theta0: f64,                 // Reference angle
    pub k_restraint: f64,
}

pub struct DihedralRestraint {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub phi0: f64,                   // Reference dihedral
    pub k_restraint: f64,
}
```

**Python Usage**: Apply biasing potentials for target geometry

---

## 6. I/O System

### 6.1 File Format Support

| Format | Type | Status | Location | Use |
|--------|------|--------|----------|-----|
| **.top** | Topology (read) | ✅ | io/topology.rs | System structure |
| **.cnf/.g96** | Coordinates (read) | ✅ | io/coordinate.rs | Initial positions |
| **.imd** | Parameters (read) | ✅ | io/imd.rs | Simulation settings |
| **.trc** | Trajectory (write) | ✅ | io/trajectory.rs | Positions/velocities/forces over time |
| **.tre** | Energy (write) | ✅ | io/energy.rs | Energy components per frame |
| **.trf** | Forces (write) | ✅ | io/force.rs | Force components per frame |
| **.ptp** | Perturbation topology (write) | ✅ | io/ptp.rs | FEP topology file |
| **.dlg** | Lambda derivatives (write) | ✅ | io/dlg.rs | dH/dλ for TI |

**Key Writing Structures**:
```rust
pub struct TrajectoryWriter {
    pub write_frame(state: &State, step: usize) -> Result<()>
}

pub struct EnergyWriter {
    pub write_frame(energies: &Energy, step: usize) -> Result<()>
}

pub struct EnergyFrame {
    pub kinetic: f64,
    pub potential: f64,
    pub total: f64,
    pub temperature: f64,
    pub pressure: f64,
    pub volume: f64,
}
```

### 6.2 Input Parsing (src/io/input/)

**GAMD Parameters**:
```rust
pub struct GamdBlock {
    pub search_mode: SearchMode,
    pub boost_form: BoostForm,
    pub threshold_type: ThresholdType,
    // ... parameters ...
}
```

**EDS Parameters**:
```rust
pub struct EdsBlock {
    pub form: EDSForm,
    pub num_states: usize,
    pub s_values: Vec<f64>,
    pub offsets: Vec<f64>,
    // ... parameters ...
}
```

**Replica Parameters**:
```rust
pub struct ReplicaBlock {
    pub num_replicas: usize,
    pub exchange_type: ExchangeType,
    pub exchange_scheme: ExchangeScheme,
    pub temperatures: Vec<f64>,
    pub lambdas: Vec<f64>,
    // ... parameters ...
}
```

---

## 7. Atom Selection System (src/selection.rs)

**Unified selection syntax** (like GROMOS++):
- `all` - All atoms
- `1-100` - Atoms 1 through 100 (1-based input, 0-based internal)
- `1,5,10-20` - Specific atoms
- `1:1-10` - Molecule 1, atoms 1-10
- `r:1-5` - Residues 1-5
- `a:CA` - All CA atoms (by name)
- `s:1-10` - Solvent molecules

**Key Class**:
```rust
pub struct AtomSelection {
    pub indices: Vec<usize>,         // 0-based
    pub n_atoms: usize,
    pub from_string(spec: &str, topo: &Topology) -> Result<Self>
    pub iter(&self) -> impl Iterator<usize>
    pub indices(&self) -> &[usize]
}
```

**Python Usage**: Select subsets for analysis or restraints

---

## 8. Analysis Tools Available (23 Command-Line Tools)

| Tool | Purpose | File |
|------|---------|------|
| **rmsd** | Root Mean Square Deviation | bin/rmsd.rs |
| **rmsf** | Root Mean Square Fluctuation | bin/rmsf.rs |
| **rgyr** | Radius of gyration | bin/rgyr.rs |
| **rdf** | Radial distribution function | bin/rdf.rs |
| **diffus** | Diffusion coefficient | bin/diffus.rs |
| **dipole** | Dipole moment analysis | bin/dipole.rs |
| **hbond** | Hydrogen bonding analysis | bin/hbond.rs |
| **ene_ana** | Energy trajectory analysis | bin/ene_ana.rs |
| **trs_ana** | Trajectory statistics | bin/trs_ana.rs |
| **frameout** | Frame extraction/formatting | bin/frameout.rs |
| **sim_box** | Box geometry analysis | bin/sim_box.rs |
| **check_top** | Topology validation | bin/check_top.rs |
| **com_top** | Center-of-mass calculations | bin/com_top.rs |
| **make_pt_top** | FEP topology generation | bin/make_pt_top.rs |
| **gamd_ana** | GaMD statistics analysis | bin/gamd_ana.rs |
| **eds_ana** | EDS analysis | bin/eds_ana.rs |
| **rep_ana** | REMD analysis | bin/rep_ana.rs |
| **pdb2g96** | PDB to GROMOS format | bin/pdb2g96.rs |

**Python Usage**: Call as subprocesses or expose core logic directly

---

## Recommended Python API Design

### Core Modules to Expose

```python
# System & State
gromos.Configuration       # System state container
gromos.State               # Single snapshot
gromos.Energy              # Energy components
gromos.Box                 # PBC box

# Structure
gromos.Topology            # Molecular system
gromos.Atom                # Individual atom

# Math
gromos.Vec3                # 3D vectors
gromos.Mat3                # 3x3 matrices

# Simulation Control
gromos.LeapFrog            # Integrator
gromos.VelocityVerlet      
gromos.StochasticDynamics
gromos.SteepestDescent

# Temperature/Pressure
gromos.BerendsenThermostat
gromos.NoséHooverThermostat
gromos.BerendsenBarostat
gromos.ParrinelloRahmanBarostat

# Constraints
gromos.SHAKE
gromos.SETTLE
gromos.LINCS

# Advanced Sampling
gromos.GamdParameters      # GaMD config
gromos.GamdRunner          # Execute GaMD
gromos.EDSParameters       # EDS config
gromos.EDSRunner           # Execute EDS
gromos.Replica             # Individual replica
gromos.ReplicaController   # REMD manager

# FEP
gromos.LambdaController    # Free energy coupling
gromos.SoftCoreParameters

# I/O
gromos.TrajectoryWriter
gromos.EnergyWriter
gromos.TopologyReader
gromos.CoordinateReader

# Utilities
gromos.AtomSelection       # Selection syntax
gromos.ForceCalculator     # Force computation
gromos.AnalysisTool        # Generic analysis framework
```

### High-Priority Exports (MVP)

1. **Configuration Management**
   - Load/save state
   - Access positions, velocities, forces
   - Energy calculation
   - Temperature/pressure computation

2. **Simulation Execution**
   - Run MD with any integrator
   - Apply thermostats/barostats
   - Handle constraints
   - Output trajectories + energies

3. **Advanced Sampling**
   - GaMD with automatic parameter search
   - EDS with multi-state sampling
   - REMD with parallel tempering
   - FEP with thermodynamic integration

4. **Force Calculations**
   - Individual bonded terms
   - Nonbonded interactions
   - Electrostatics (CRF, PME)
   - Restraint potentials

5. **I/O Operations**
   - Read topology + coordinates
   - Write trajectories + energies
   - FEP file handling
   - Parameter file parsing

### Medium-Priority Exports

6. **Analysis Tools**
   - RMSD, RMSF, RDF calculations
   - Trajectory statistics
   - Energy analysis
   - Post-processing utilities

7. **Detailed Control**
   - Individual force components
   - Energy group accounting
   - Pairlist management
   - Selection system

### Nice-to-Have

8. **Direct C API**
   - FFI for existing C++ code
   - In-process library linking

---

## Performance Characteristics

### Optimizations Already Implemented

1. **SIMD Vectorization**
   - glam (v0.24) with bytemuck for Vec3/Mat3
   - wide crate for wide SIMD operations
   - Automatic compiler vectorization

2. **Parallelization**
   - Rayon for force calculation parallelism
   - Crossbeam for thread management
   - Lock-free where possible

3. **Memory Efficiency**
   - MiMalloc allocator (optional, feature flag)
   - Double-buffering for state (O(1) pointer swap)
   - Grid-cell pairlist (O(N) vs O(N²))

4. **Algorithm Choices**
   - Leap-Frog integrator (good energy conservation)
   - Reaction Field (GROMOS standard, fast)
   - Cell-neighbor approach (vs full pairlist)

### Benchmarking Features

Built-in diagnostics:
- EnergyDiagnostics
- ForceDiagnostics  
- FrameDiagnostics
- ProgressBar
- Timer utilities

---

## Testing Infrastructure

| Type | Files | Coverage |
|------|-------|----------|
| Unit tests | In each .rs | Force calculations, energy |
| Integration tests | tests/md_integration_tests.rs | Full MD workflow |
| Advanced tests | tests/md_advanced_features_tests.rs | GaMD, EDS, FEP |
| Binary tests | tests/md_binary_tests.rs | CLI interface |

---

## Implementation Roadmap for Python Binding

### Phase 1: Core Classes (1-2 weeks)
- [ ] Vec3, Mat3 - basic math types
- [ ] Configuration, State, Energy - system state
- [ ] Topology, Atom - structure definition
- [ ] Load/save topology, coordinates

### Phase 2: Basic Simulation (1-2 weeks)
- [ ] Integrators (LeapFrog, VelocityVerlet, SD)
- [ ] Thermostats + Barostats
- [ ] Force calculation (bonded + nonbonded)
- [ ] Basic MD loop

### Phase 3: Constraints & Advanced (2-3 weeks)
- [ ] SHAKE, SETTLE, LINCS
- [ ] Reaction Field electrostatics
- [ ] Restraints
- [ ] Output (trajectory, energy, forces)

### Phase 4: Advanced Sampling (3-4 weeks)
- [ ] FEP with soft-core
- [ ] GAMD implementation
- [ ] EDS implementation
- [ ] REMD implementation

### Phase 5: Analysis Tools (2-3 weeks)
- [ ] Expose calculation functions (RMSD, RMSF, RDF)
- [ ] Higher-level analysis objects
- [ ] Post-processing utilities

### Phase 6: Polish (1-2 weeks)
- [ ] Documentation + examples
- [ ] Error handling
- [ ] Performance optimization
- [ ] PyPI packaging

---

## Key Design Decisions

### 1. What to Expose vs Hide
**EXPOSE**: High-level simulation control, state management, I/O
**HIDE**: Internal pairlist structures, memory management details, threading primitives

### 2. Data Conversion Strategy
**Efficient**: Use NumPy's array protocol where possible
**Safe**: Implement bounds checking, type validation
**Fast**: Minimize copies between Rust and Python

### 3. Error Handling
- All Result types → Python exceptions
- Validation at module boundaries
- Clear error messages

### 4. Thread Safety
- Configuration/State structures wrapped in Arc<Mutex<...>> if needed
- Python GIL considerations
- Parallel force calculations via rayon (no GIL release needed if careful)

---

## Conclusion

The **gromos-rs** codebase provides a comprehensive, production-ready MD engine with advanced sampling capabilities. A Python API should prioritize:

1. **Core simulation workflow** (load → setup → run → analyze)
2. **Easy-to-use integrators** and **thermostats**
3. **Advanced sampling methods** (GAMD, EDS, REMD)
4. **Flexible I/O** and **analysis tools**
5. **Good error messages** and **examples**

The modular Rust design makes this achievable with PyO3 bindings, with clear module boundaries and minimal abstraction overhead.

