# GROMOS-RS Python API Reference

Complete API documentation for the GROMOS-RS Python bindings.

## Table of Contents

1. [Math Types](#math-types)
   - [Vec3](#vec3)
   - [Mat3](#mat3)
2. [Core Structures](#core-structures)
   - [Box](#box)
   - [Energy](#energy)
   - [State](#state)
   - [Configuration](#configuration)
   - [Topology](#topology)
3. [Integrators](#integrators)
   - [LeapFrog](#leapfrog)
   - [VelocityVerlet](#velocityverlet)
   - [StochasticDynamics](#stochasticdynamics)
4. [Advanced Sampling](#advanced-sampling)
   - [GaMD](#gamd)
   - [EDS](#eds)
   - [REMD](#remd)

---

## Math Types

### Vec3

**SIMD-accelerated 3D vector**

Zero-copy conversion to/from NumPy arrays. All operations use SIMD instructions for performance.

#### Constructor

```python
Vec3(x: float, y: float, z: float) -> Vec3
```

Create a new 3D vector.

**Parameters:**
- `x` (float): X component
- `y` (float): Y component
- `z` (float): Z component

**Example:**
```python
v = gromos.Vec3(1.0, 2.0, 3.0)
```

#### Properties

- `x: float` - X component (read/write)
- `y: float` - Y component (read/write)
- `z: float` - Z component (read/write)

**Example:**
```python
v = gromos.Vec3(1.0, 2.0, 3.0)
print(v.x)  # 1.0
v.x = 5.0
```

#### Methods

##### `length() -> float`

Vector magnitude (L2 norm).

**Example:**
```python
v = gromos.Vec3(3.0, 4.0, 0.0)
print(v.length())  # 5.0
```

##### `length_squared() -> float`

Squared length (faster than `length()`).

**Example:**
```python
v = gromos.Vec3(3.0, 4.0, 0.0)
print(v.length_squared())  # 25.0
```

##### `normalize() -> Vec3`

Return normalized (unit length) vector.

**Example:**
```python
v = gromos.Vec3(3.0, 4.0, 0.0)
v_norm = v.normalize()
print(v_norm.length())  # 1.0
```

##### `dot(other: Vec3) -> float`

Dot product with another vector.

**Example:**
```python
v1 = gromos.Vec3(1.0, 0.0, 0.0)
v2 = gromos.Vec3(0.0, 1.0, 0.0)
print(v1.dot(v2))  # 0.0
```

##### `cross(other: Vec3) -> Vec3`

Cross product with another vector.

**Example:**
```python
v1 = gromos.Vec3(1.0, 0.0, 0.0)
v2 = gromos.Vec3(0.0, 1.0, 0.0)
v3 = v1.cross(v2)  # Vec3(0.0, 0.0, 1.0)
```

##### `distance(other: Vec3) -> float`

Euclidean distance to another vector.

**Example:**
```python
v1 = gromos.Vec3(0.0, 0.0, 0.0)
v2 = gromos.Vec3(3.0, 4.0, 0.0)
print(v1.distance(v2))  # 5.0
```

##### `to_numpy() -> np.ndarray`

Convert to NumPy array (shape: (3,), dtype: float32).

**Example:**
```python
v = gromos.Vec3(1.0, 2.0, 3.0)
arr = v.to_numpy()
print(arr)  # [1. 2. 3.]
```

##### `from_numpy(arr: np.ndarray) -> Vec3` (static)

Create from NumPy array.

**Parameters:**
- `arr` (np.ndarray): Array with 3 elements

**Example:**
```python
arr = np.array([1.0, 2.0, 3.0], dtype=np.float32)
v = gromos.Vec3.from_numpy(arr)
```

#### Operators

- `v1 + v2` - Vector addition
- `v1 - v2` - Vector subtraction
- `v * scalar` - Scalar multiplication

**Example:**
```python
v1 = gromos.Vec3(1.0, 2.0, 3.0)
v2 = gromos.Vec3(4.0, 5.0, 6.0)
v3 = v1 + v2
v4 = v1 * 2.0
```

---

### Mat3

**SIMD-accelerated 3×3 matrix**

#### Constructor

##### `identity() -> Mat3` (static)

Create identity matrix.

**Example:**
```python
mat = gromos.Mat3.identity()
```

##### `from_cols(x: Vec3, y: Vec3, z: Vec3) -> Mat3` (static)

Create matrix from column vectors.

**Example:**
```python
v1 = gromos.Vec3(1.0, 0.0, 0.0)
v2 = gromos.Vec3(0.0, 1.0, 0.0)
v3 = gromos.Vec3(0.0, 0.0, 1.0)
mat = gromos.Mat3.from_cols(v1, v2, v3)
```

#### Methods

##### `determinant() -> float`

Matrix determinant.

##### `inverse() -> Mat3`

Matrix inverse.

##### `transpose() -> Mat3`

Matrix transpose.

##### `mul_vec3(v: Vec3) -> Vec3`

Matrix-vector multiplication.

**Example:**
```python
mat = gromos.Mat3.identity()
v = gromos.Vec3(1.0, 2.0, 3.0)
result = mat.mul_vec3(v)
```

##### `to_numpy() -> np.ndarray`

Convert to NumPy array (shape: (3, 3), dtype: float32).

---

## Core Structures

### Box

**Simulation box with periodic boundaries**

#### Constructors

##### `vacuum() -> Box` (static)

Create vacuum box (no periodicity).

**Example:**
```python
box = gromos.Box.vacuum()
```

##### `rectangular(lx: float, ly: float, lz: float) -> Box` (static)

Create rectangular box with dimensions in nm.

**Parameters:**
- `lx` (float): Length along x-axis (nm)
- `ly` (float): Length along y-axis (nm)
- `lz` (float): Length along z-axis (nm)

**Example:**
```python
box = gromos.Box.rectangular(3.0, 3.0, 3.0)  # 3×3×3 nm box
```

##### `triclinic(mat: Mat3) -> Box` (static)

Create triclinic box from box vectors.

**Parameters:**
- `mat` (Mat3): Box vectors as columns

**Example:**
```python
# Create triclinic box
v1 = gromos.Vec3(3.0, 0.0, 0.0)
v2 = gromos.Vec3(0.5, 3.0, 0.0)
v3 = gromos.Vec3(0.0, 0.0, 3.0)
mat = gromos.Mat3.from_cols(v1, v2, v3)
box = gromos.Box.triclinic(mat)
```

#### Methods

##### `volume() -> float`

Box volume in nm³.

**Example:**
```python
box = gromos.Box.rectangular(3.0, 3.0, 3.0)
print(box.volume())  # 27.0
```

##### `dimensions() -> Vec3`

Box dimensions (for rectangular box).

**Example:**
```python
box = gromos.Box.rectangular(3.0, 4.0, 5.0)
dims = box.dimensions()
print(dims.x, dims.y, dims.z)  # 3.0 4.0 5.0
```

---

### Energy

**Energy storage for molecular system**

All energies in kJ/mol.

#### Constructor

```python
Energy(num_temperature_groups: int, num_energy_groups: int) -> Energy
```

Create energy object.

**Parameters:**
- `num_temperature_groups` (int): Number of temperature groups
- `num_energy_groups` (int): Number of energy groups

**Example:**
```python
energy = gromos.Energy(num_temperature_groups=1, num_energy_groups=1)
```

#### Properties (read-only)

- `kinetic: float` - Total kinetic energy
- `potential: float` - Total potential energy
- `bond: float` - Bond stretch energy
- `angle: float` - Angle bend energy
- `dihedral: float` - Dihedral torsion energy
- `lj: float` - Lennard-Jones energy
- `coulomb: float` - Coulomb (electrostatic) energy

#### Methods

##### `total() -> float`

Total energy (kinetic + potential).

**Example:**
```python
print(energy.total())
```

##### `clear()`

Clear all energies to zero.

**Example:**
```python
energy.clear()
```

##### `to_dict() -> dict`

Get all energies as dictionary.

**Example:**
```python
d = energy.to_dict()
print(d['kinetic'])
print(d['potential'])
```

---

### State

**System state (positions, velocities, forces)**

Zero-copy data sharing with NumPy arrays.

#### Constructor

```python
State(num_atoms: int, num_temp_groups: int, num_energy_groups: int) -> State
```

Create state for N atoms.

**Parameters:**
- `num_atoms` (int): Number of atoms
- `num_temp_groups` (int): Number of temperature groups
- `num_energy_groups` (int): Number of energy groups

**Example:**
```python
state = gromos.State(
    num_atoms=1000,
    num_temp_groups=1,
    num_energy_groups=1
)
```

#### Methods

##### `num_atoms() -> int`

Number of atoms.

##### `positions() -> np.ndarray`

Get positions as NumPy array (N × 3, float32).

**Returns:** Array of shape (N, 3) in nm

**Example:**
```python
pos = state.positions()
print(pos.shape)  # (1000, 3)
```

##### `velocities() -> np.ndarray`

Get velocities as NumPy array (N × 3, float32).

**Returns:** Array of shape (N, 3) in nm/ps

##### `forces() -> np.ndarray`

Get forces as NumPy array (N × 3, float32).

**Returns:** Array of shape (N, 3) in kJ/(mol·nm)

##### `set_positions(arr: np.ndarray)`

Set positions from NumPy array.

**Parameters:**
- `arr` (np.ndarray): Array of shape (N, 3), dtype float32

**Example:**
```python
pos = np.random.rand(1000, 3).astype(np.float32)
state.set_positions(pos)
```

##### `set_velocities(arr: np.ndarray)`

Set velocities from NumPy array.

**Parameters:**
- `arr` (np.ndarray): Array of shape (N, 3), dtype float32

---

### Configuration

**Complete system configuration**

Combines state, energy, and topology.

#### Constructor

```python
Configuration(num_atoms: int, num_temp_groups: int, num_energy_groups: int) -> Configuration
```

#### Methods

##### `current_state() -> State`

Get current state.

##### `current_energy() -> Energy`

Get current energy.

---

### Topology

**Molecular topology (atoms, bonds, parameters)**

#### Constructor

```python
Topology() -> Topology
```

Create empty topology.

**Example:**
```python
topo = gromos.Topology()
```

#### Methods

##### `num_atoms() -> int`

Number of atoms.

##### `num_bonds() -> int`

Number of bonds.

##### `num_angles() -> int`

Number of angles.

##### `num_dihedrals() -> int`

Number of dihedrals.

---

## Integrators

### LeapFrog

**Leap-Frog integrator (velocity Verlet variant)**

Fast and stable. Velocities offset by dt/2 from positions.

#### Constructor

```python
LeapFrog(dt: float) -> LeapFrog
```

**Parameters:**
- `dt` (float): Timestep in ps

**Example:**
```python
integrator = gromos.LeapFrog(dt=0.002)  # 2 fs
```

#### Methods

##### `timestep() -> float`

Get timestep in ps.

---

### VelocityVerlet

**Velocity Verlet integrator**

Higher accuracy than Leap-Frog. Positions and velocities at same time.

#### Constructor

```python
VelocityVerlet(dt: float) -> VelocityVerlet
```

**Parameters:**
- `dt` (float): Timestep in ps

**Example:**
```python
integrator = gromos.VelocityVerlet(dt=0.001)  # 1 fs
```

#### Methods

##### `timestep() -> float`

Get timestep in ps.

---

### StochasticDynamics

**Stochastic Dynamics (Langevin) integrator**

Implicit solvent with friction and random forces.

#### Constructor

```python
StochasticDynamics(dt: float, gamma: float, temperature: float) -> StochasticDynamics
```

**Parameters:**
- `dt` (float): Timestep in ps
- `gamma` (float): Friction coefficient in ps⁻¹
- `temperature` (float): Target temperature in K

**Example:**
```python
integrator = gromos.StochasticDynamics(
    dt=0.002,
    gamma=0.1,
    temperature=300.0
)
```

#### Methods

##### `timestep() -> float`

Get timestep in ps.

---

## Advanced Sampling

### GaMD

**Gaussian Accelerated Molecular Dynamics**

Adds harmonic boost potential to smooth energy landscape.

#### GamdParameters

```python
GamdParameters(sigma0: float, threshold_mode: str) -> GamdParameters
```

**Parameters:**
- `sigma0` (float): Standard deviation of boost (kJ/mol)
- `threshold_mode` (str): 'lower' or 'upper'

**Example:**
```python
params = gromos.GamdParameters(sigma0=6.0, threshold_mode='lower')
```

#### GamdRunner

```python
GamdRunner(params: GamdParameters) -> GamdRunner
```

**Example:**
```python
runner = gromos.GamdRunner(params)
```

---

### EDS

**Enveloping Distribution Sampling**

Samples multiple end-states simultaneously with smoothed envelope.

#### EDSParameters

```python
EDSParameters(num_states: int, smoothness: float) -> EDSParameters
```

**Parameters:**
- `num_states` (int): Number of end-states
- `smoothness` (float): Envelope smoothness (kJ/mol)

**Example:**
```python
params = gromos.EDSParameters(num_states=4, smoothness=1.0)
```

#### EDSRunner

```python
EDSRunner(params: EDSParameters) -> EDSRunner
```

**Example:**
```python
runner = gromos.EDSRunner(params)
```

---

### REMD

**Replica Exchange Molecular Dynamics**

Manages multiple replicas with periodic exchange attempts.

#### ReplicaController

```python
ReplicaController(num_replicas: int, exchange_interval: int) -> ReplicaController
```

**Parameters:**
- `num_replicas` (int): Number of replicas
- `exchange_interval` (int): Steps between exchange attempts

**Example:**
```python
remd = gromos.ReplicaController(
    num_replicas=8,
    exchange_interval=1000
)
```

#### Methods

##### `num_replicas() -> int`

Get number of replicas.

---

## Units

GROMOS-RS uses the following units:

| Quantity | Unit |
|----------|------|
| Length | nm (nanometers) |
| Time | ps (picoseconds) |
| Energy | kJ/mol |
| Force | kJ/(mol·nm) |
| Velocity | nm/ps |
| Temperature | K (Kelvin) |
| Mass | g/mol (atomic mass units) |

## Performance Tips

1. **Zero-copy NumPy**: Use `.positions()`, `.velocities()`, `.forces()` to get views
2. **SIMD acceleration**: Vec3 operations automatically vectorized
3. **Parallel execution**: Force calculations parallelized via Rayon
4. **Batch operations**: Process multiple atoms at once when possible

## Examples

See the `examples/` directory for complete examples:

- `01_basic_vectors.py` - Vector and matrix operations
- `02_system_setup.py` - System initialization
- `03_integrators.py` - Integration algorithms
- `04_advanced_sampling.py` - Enhanced sampling methods
- `05_complete_workflow.py` - Full simulation workflow

## See Also

- [README.md](README.md) - Installation and quick start
- [GROMOS documentation](http://www.gromos.net)
- [PyO3 documentation](https://pyo3.rs)
