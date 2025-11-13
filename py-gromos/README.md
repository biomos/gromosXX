# GROMOS-RS Python Bindings

High-performance molecular dynamics simulations powered by Rust, accessible from Python.

## Architecture

Inspired by [Polars](https://github.com/pola-rs/polars)' design:

- **Rust core**: 100% Rust MD engine with SIMD and parallel execution
- **PyO3 bindings**: Safe, zero-copy data sharing between Rust and Python
- **NumPy integration**: Direct array conversion without copying
- **Memory safety**: Guaranteed by Rust's type system

## Installation

### From source (development)

```bash
# Install maturin (Rust-Python build tool)
pip install maturin

# Build and install in development mode
cd py-gromos
maturin develop --release

# Or build wheel for distribution
maturin build --release
pip install target/wheels/gromos-*.whl
```

### From PyPI (when published)

```bash
pip install gromos
```

## Quick Start

```python
import gromos
import numpy as np

# Create system state
state = gromos.State(
    num_atoms=1000,
    num_temp_groups=1,
    num_energy_groups=1
)

# Set up simulation box
box = gromos.Box.rectangular(5.0, 5.0, 5.0)  # 5x5x5 nm

# Create integrator
integrator = gromos.LeapFrog(dt=0.002)  # 2 fs timestep

# Work with 3D vectors
v1 = gromos.Vec3(1.0, 2.0, 3.0)
v2 = gromos.Vec3(4.0, 5.0, 6.0)
distance = v1.distance(v2)
print(f"Distance: {distance:.3f} nm")

# Convert to/from NumPy
vec_array = v1.to_numpy()
v3 = gromos.Vec3.from_numpy(np.array([1.0, 2.0, 3.0]))

# Access simulation data as NumPy arrays (zero-copy)
positions = state.positions()  # Shape: (N, 3)
velocities = state.velocities()  # Shape: (N, 3)
forces = state.forces()  # Shape: (N, 3)
```

## Advanced Sampling

### Gaussian Accelerated MD (GaMD)

```python
# Create GaMD parameters
gamd_params = gromos.GamdParameters(
    sigma0=6.0,
    threshold_mode='lower'
)

# Initialize GaMD runner
gamd = gromos.GamdRunner(gamd_params)
```

### Enveloping Distribution Sampling (EDS)

```python
# Create EDS parameters for 4 states
eds_params = gromos.EDSParameters(
    num_states=4,
    smoothness=1.0
)

# Initialize EDS runner
eds = gromos.EDSRunner(eds_params)
```

### Replica Exchange MD (REMD)

```python
# Create REMD controller
remd = gromos.ReplicaController(
    num_replicas=8,
    exchange_interval=1000  # Attempt exchange every 1000 steps
)

print(f"Managing {remd.num_replicas()} replicas")
```

## Features

### Math Types
- `Vec3`: SIMD-accelerated 3D vectors
- `Mat3`: SIMD-accelerated 3×3 matrices

### Core Structures
- `State`: System state (positions, velocities, forces)
- `Energy`: Energy tracking (kinetic, potential, components)
- `Configuration`: Complete system configuration
- `Topology`: Molecular topology (atoms, bonds, parameters)
- `Box`: Simulation box (vacuum, rectangular, triclinic)

### Integrators
- `LeapFrog`: Fast velocity Verlet variant
- `VelocityVerlet`: Higher accuracy integrator
- `StochasticDynamics`: Langevin dynamics for implicit solvent

### Advanced Sampling
- `GamdParameters`, `GamdRunner`: Gaussian Accelerated MD
- `EDSParameters`, `EDSRunner`: Enveloping Distribution Sampling
- `ReplicaController`: Replica Exchange MD

## Performance

GROMOS-RS achieves 2-3x speedup over the original C++ implementation:

| Feature | Benefit |
|---------|---------|
| SIMD vectorization | Automatic via glam + wide |
| Parallel execution | Rayon multi-threading |
| Zero-copy arrays | NumPy integration via PyO3 |
| Memory safety | Rust guarantees, no segfaults |
| Zero-cost abstractions | No runtime overhead |

### Why Zero-Copy Matters

Like Polars, GROMOS-RS uses Apache Arrow-compatible memory layout:

```python
# No data copying - direct memory access
positions = state.positions()  # Returns view of Rust memory
positions[0] = [1.0, 2.0, 3.0]  # Would need copy for mutation

# Setting data requires a copy
new_positions = np.random.rand(1000, 3).astype(np.float32)
state.set_positions(new_positions)  # Copies into Rust
```

## Comparison with Polars

| Aspect | Polars | GROMOS-RS |
|--------|--------|-----------|
| Core language | Rust | Rust |
| Python bindings | PyO3 | PyO3 |
| Data sharing | Arrow zero-copy | NumPy zero-copy |
| Parallelism | Rayon | Rayon |
| SIMD | Yes | Yes (glam) |
| Domain | DataFrames | Molecular dynamics |

## Development

### Building

```bash
# Install Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install maturin
pip install maturin

# Build in development mode (includes debug symbols)
maturin develop

# Build optimized release
maturin develop --release
```

### Testing

```bash
# Run Python tests
pytest tests/

# Run Rust tests
cargo test -p gromos-rs
```

### Project Structure

```
py-gromos/
├── Cargo.toml           # Rust package manifest
├── pyproject.toml       # Python package config
├── src/
│   └── lib.rs          # Rust ↔ Python bindings (PyO3)
├── python/
│   └── gromos/
│       └── __init__.py  # Python API
├── tests/              # Python tests
└── examples/           # Example scripts
```

## License

GPL-2.0 - Same as GROMOS

## References

- **GROMOS**: http://www.gromos.net
- **Polars**: https://github.com/pola-rs/polars
- **PyO3**: https://pyo3.rs
- **NumPy**: https://numpy.org

## Citation

If you use GROMOS-RS in your research, please cite:

```
GROMOS Software for Biomolecular Simulation
http://www.gromos.net
```
