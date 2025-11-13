# Python API for gromos-rs: Executive Summary

## Overview

**gromos-rs** is a high-performance Rust implementation of GROMOS molecular dynamics with:
- **7,500+ lines** of optimized Rust code
- **Production-ready core MD** engine
- **Advanced sampling**: GAMD, EDS, REMD
- **23 analysis tools** for trajectory post-processing
- **SIMD + parallelization** for performance

---

## What Should Be Exposed to Python

### TIER 1: Core Functionality (Must-Have)

**1. System State Management**
```python
# Load and manage molecular system
topology = gromos.Topology.from_file("system.top")
configuration = gromos.Configuration(n_atoms=1000)

# Access state snapshots
positions = configuration.current().pos        # List[Vec3]
velocities = configuration.current().vel
forces = configuration.current().force

# Calculate properties
kinetic_energy = configuration.current().calculate_kinetic_energy(masses)
temperature = configuration.current().temperature(n_dof)
pressure = configuration.current().pressure()
```

**2. Basic Molecular Dynamics**
```python
# Choose integrator
integrator = gromos.LeapFrog()  # or VelocityVerlet, StochasticDynamics

# Control temperature & pressure
thermostat = gromos.BerendsenThermostat(temperature=300.0, tau_t=0.1)
barostat = gromos.BerendsenBarostat(pressure=1.0, tau_p=0.5)

# Apply constraints
constraints = gromos.SHAKE(tolerance=1e-4)  # or SETTLE, LINCS

# Run simulation
for step in range(1000):
    integrator.step(dt=0.002, topology=topology, configuration=configuration)
    thermostat.apply(configuration=configuration, masses=masses)
    forces = calculate_forces(topology, configuration)
```

**3. I/O Operations**
```python
# Read inputs
topo = gromos.read_topology("system.top")
coords = gromos.read_coordinates("initial.g96")

# Write outputs
traj_writer = gromos.TrajectoryWriter("output.trc")
energy_writer = gromos.EnergyWriter("output.tre")

# Write trajectory frames
traj_writer.write_frame(configuration.current(), step=0)
energy_writer.write_frame(configuration.current().energies, step=0)
```

---

### TIER 2: Advanced Sampling (Important)

**1. Gaussian Accelerated MD (GaMD)**
```python
# Automatic parameter adaptation
gamd_params = gromos.GamdParameters(
    search_mode=gromos.SearchMode.CmdSearch,  # Data collection
    boost_form=gromos.BoostForm.DualBoost,
    threshold_type=gromos.ThresholdType.LowerBound,
    sigma0_dih=6.0,
    sigma0_tot=12.0,
)

# Run simulation with automatic statistics
runner = gromos.GamdRunner(gamd_params, topology, configuration)
runner.run(n_steps=10000, dt=0.002, integrator=integrator)

# Access statistics
stats = runner.statistics()
print(f"Dihedral V_max={stats.dih_v_max}, V_mean={stats.dih_v_mean}")
```

**2. Enveloping Distribution Sampling (EDS)**
```python
# Multi-state sampling
eds_params = gromos.EDSParameters(
    form=gromos.EDSForm.SingleS,
    num_states=2,
    s_values=[0.5],
    offsets=[0.0, 1000.0],  # State B offset
    temperature=300.0,
    n_atoms=topology.num_atoms(),
)

# Run EDS
runner = gromos.EDSRunner(eds_params, topology, configuration)
runner.run(n_steps=10000, dt=0.002)

# Query sampling statistics
current_state = runner.current_state_id()
round_trips = runner.round_trips()
```

**3. Replica Exchange MD (REMD)**
```python
# Setup replicas at different temperatures
temperatures = [300, 320, 340, 360]
replicas = [
    gromos.Replica(temperature=T) 
    for T in temperatures
]

# Control exchange
controller = gromos.ReplicaController(
    replicas=replicas,
    exchange_type=gromos.ExchangeType.Temperature,
    exchange_scheme=gromos.ExchangeScheme.OddEven,
    exchange_interval=100,
)

# Run parallel replicas
controller.run_parallel(n_steps=10000, dt=0.002)

# Check acceptance statistics
stats = controller.statistics()
for i in range(len(temperatures)):
    for j in range(i+1, len(temperatures)):
        rate = stats.pair_acceptance_rate(i, j)
        print(f"T{i}-T{j}: {rate:.2%} acceptance")
```

---

### TIER 3: Free Energy Perturbation (Nice-to-Have)

**1. Thermodynamic Integration**
```python
# Setup FEP
fep = gromos.LambdaController()
fep.set_lambda(0.0)  # Start in state A
fep.set_dlambda_dt(0.00001)  # Slow growth

# Run slow-growth FEP
for step in range(100000):
    fep.update_lambda(dt=0.002)
    
    # Calculate forces with current lambda
    forces_A = calculate_forces_state_A(...)
    forces_B = calculate_forces_state_B(...)
    forces = interpolate_forces(forces_A, forces_B, fep.lambda)
    
    # Accumulate dH/dlambda
    dH_dlambda = compute_free_energy_derivative(forces_A, forces_B, fep.lambda)
```

**2. Soft-Core Potentials**
```python
# Prevent singularities in FEP
softcore = gromos.SoftCoreParameters(
    alpha=0.5,      # Soft-core exponent
    lambda_power=1,
    sc_sigma=0.3,   # Min distance in nm
)
```

---

## Key Data Structures

### Math Types
```python
vec = gromos.Vec3(x=1.0, y=2.0, z=3.0)
mat = gromos.Mat3(...)  # 3x3 matrix

# Operations
dot_product = vec1.dot(vec2)
cross_product = vec1.cross(vec2)
length = vec.length()
```

### Molecular Structure
```python
atom = topology.atoms[0]  # gromos.Atom
print(f"{atom.name} (mass={atom.mass}, charge={atom.charge})")

# Bonds, angles, dihedrals
for bond in topology.bonds:
    print(f"Bond {bond.i}-{bond.j}: type {bond.bond_type}")
```

### Energy Components
```python
energy = configuration.current().energies

print(f"Kinetic:    {energy.kinetic_total:.2f} kJ/mol")
print(f"Bond:       {energy.bond_total:.2f} kJ/mol")
print(f"Angle:      {energy.angle_total:.2f} kJ/mol")
print(f"Dihedral:   {energy.dihedral_total:.2f} kJ/mol")
print(f"LJ:         {energy.lj_total:.2f} kJ/mol")
print(f"Coulomb:    {energy.crf_total:.2f} kJ/mol")
print(f"Total:      {energy.total():.2f} kJ/mol")
```

### Atom Selection
```python
# Unified selection syntax (like GROMOS++)
selection = gromos.AtomSelection("a:CA", topology)      # All CA atoms
selection = gromos.AtomSelection("r:1-10", topology)    # Residues 1-10
selection = gromos.AtomSelection("1:1-50", topology)    # Mol 1, atoms 1-50

# Use in analysis
for idx in selection:
    print(f"Selected atom {idx}")
```

---

## Design Principles

### 1. Python First
- Pythonic naming (snake_case for functions/methods)
- NumPy compatibility (array conversions)
- Intuitive class hierarchy

### 2. Zero-Copy Where Possible
- Pass references, not copies
- Use buffer protocol for arrays
- Lazy evaluation where appropriate

### 3. Sensible Defaults
```python
# Simple case
md = gromos.MD(topology, configuration)
md.run(n_steps=1000)  # Uses defaults

# Advanced case
md.set_integrator(gromos.VelocityVerlet())
md.set_thermostat(gromos.NoséHooverThermostat(300.0))
md.set_constraints(gromos.SHAKE())
```

### 4. Rich Error Messages
```python
try:
    config = gromos.Configuration(n_atoms=-1)
except ValueError as e:
    print(e)  # "n_atoms must be > 0, got -1"
```

---

## Priority Implementation Order

### Phase 1 (MVP - 1-2 weeks)
1. Math types (Vec3, Mat3)
2. Configuration, State, Energy
3. Topology loading
4. Basic I/O (trajectory, energy)

### Phase 2 (Basic MD - 1-2 weeks)
1. Integrators
2. Thermostats + barostats
3. Force calculations (bonded + nonbonded)
4. Constraints (SHAKE)

### Phase 3 (Advanced - 2-3 weeks)
1. GAMD
2. EDS
3. REMD
4. Additional constraints (SETTLE, LINCS)

### Phase 4 (Polish - 1-2 weeks)
1. Documentation
2. Examples
3. Performance optimization
4. PyPI packaging

---

## File Locations in Codebase

| Component | File(s) | Lines |
|-----------|---------|-------|
| Configuration/Energy | src/configuration.rs | 484 |
| Topology | src/topology.rs | 1,297 |
| Math (Vec3, Mat3) | Re-exported from glam | - |
| Integrators | src/integrator.rs | 694 |
| Thermostats/Barostats | src/algorithm/ | 350+ |
| Constraints | src/algorithm/constraints.rs | 770 |
| Bonded Forces | src/interaction/bonded.rs | 2,624 |
| Nonbonded Forces | src/interaction/nonbonded.rs | 907 |
| Electrostatics | src/interaction/electrostatics.rs | 1,137 |
| Restraints | src/interaction/restraints.rs | 858 |
| GAMD | src/gamd.rs | 864 |
| EDS | src/eds.rs | 1,036 |
| REMD | src/remd.rs + src/replica.rs | 1,277 + ? |
| FEP | src/fep.rs | 687 |
| I/O | src/io/ | ~3,500 |

---

## Example: Complete Simulation Workflow

```python
import gromos

# 1. Load system
topology = gromos.Topology.from_file("system.top")
config = gromos.Configuration.from_file("initial.g96")

# 2. Setup simulation parameters
integrator = gromos.LeapFrog()
thermostat = gromos.BerendsenThermostat(temperature=300.0, tau_t=0.1)
constraints = gromos.SHAKE(tolerance=1e-4)

# 3. Setup output
traj_writer = gromos.TrajectoryWriter("output.trc")
energy_writer = gromos.EnergyWriter("output.tre")

# 4. Run simulation
for step in range(10000):
    # Calculate forces
    bonded_forces = gromos.calculate_bonded_forces(topology, config)
    nonbonded_forces = gromos.calculate_nonbonded_forces(topology, config)
    
    # Integrate
    integrator.step(dt=0.002, topology=topology, config=config)
    
    # Thermostat
    thermostat.apply(config)
    
    # Constraints
    constraints.apply(topology, config)
    
    # Output every 100 steps
    if step % 100 == 0:
        traj_writer.write_frame(config.current(), step)
        energy_writer.write_frame(config.current().energies, step)
        
        T = config.current().temperature(3*topology.num_atoms() - len(constraints.bonds))
        E = config.current().energies.total()
        print(f"Step {step}: T={T:.1f} K, E={E:.2f} kJ/mol")

print("Simulation complete!")
```

---

## Technical Considerations

### Python-Rust Bridge Strategy: PyO3

**Advantages**:
- Native Rust semantics
- Minimal overhead
- GIL-release capable
- NumPy integration

**Module structure**:
```
gromos_py/              # Python package
├── __init__.py         # Main entry point
├── _gromos.so          # Compiled extension
└── stubs/              # Type hints (.pyi)
```

### Memory Management
- Rust handles all memory (no manual Python allocation)
- Reference counting via Arc<> where shared
- Configuration/State can be large (1000s of atoms × 3 coordinates × 8 bytes)

### Serialization
- Use GROMOS file formats (.top, .cnf, .trc, .tre)
- Optional NumPy array export for Python analysis

---

## Estimated Effort

| Component | Complexity | Time | Status |
|-----------|-----------|------|--------|
| Math types | Low | 1-2 days | Foundation |
| Core data structures | Low | 2-3 days | Foundation |
| Basic I/O | Low | 2-3 days | Foundation |
| Integrators | Medium | 3-4 days | Integration |
| Forces | High | 1-2 weeks | Integration |
| GAMD | High | 1 week | Advanced |
| EDS | High | 1 week | Advanced |
| REMD | High | 1-2 weeks | Advanced |
| Analysis tools | Medium | 1 week | Post-processing |
| Documentation | Medium | 1 week | Polish |
| **TOTAL** | - | **6-8 weeks** | |

---

## References

- **Feature completeness**: See `/home/user/gromosXX/GROMOS_FEATURES_CATALOG.md`
- **Full technical analysis**: See `/home/user/gromosXX/PYTHON_API_DESIGN.md`
- **Source code**: `/home/user/gromosXX/gromos-rs/src/`

