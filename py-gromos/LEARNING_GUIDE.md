# GROMOS-RS Python Bindings - Learning Guide

Comprehensive learning materials to understand the Polars-inspired architecture.

## 📚 Overview

This guide provides multiple paths to learn about GROMOS-RS Python bindings:
- **Notebooks**: Interactive Jupyter notebooks (theory + practice)
- **Scripts**: Standalone Python scripts (hands-on demos)
- **Examples**: Production-ready code samples

## 🎯 Learning Paths

### Path 1: Quick Start (30 minutes)

For users who want to start using the bindings quickly:

1. **Read**: `README.md` - Installation and basics
2. **Run**: `scripts/hands_on_demo.py` - Interactive demos
3. **Try**: `examples/01_basic_vectors.py` - Simple example

```bash
cd py-gromos
python scripts/hands_on_demo.py
python examples/01_basic_vectors.py
```

### Path 2: Deep Understanding (2-3 hours)

For users who want to understand the architecture:

1. **Read**: `scripts/explain_architecture.py` - Architecture explanation
2. **Study**: `notebooks/01_understanding_pyo3_bindings.ipynb` - PyO3 concepts
3. **Explore**: `notebooks/02_molecular_systems_and_energy.ipynb` - MD systems
4. **Optimize**: `notebooks/03_performance_deep_dive.ipynb` - Performance

```bash
cd py-gromos
python scripts/explain_architecture.py
jupyter notebook notebooks/
```

### Path 3: Complete Mastery (1 day)

For developers who want to extend or contribute:

1. All materials from Path 2
2. Read the source: `src/lib.rs` - Rust implementation
3. Study examples: All files in `examples/`
4. API reference: `API_REFERENCE.md`

## 📖 Materials Index

### Notebooks (`notebooks/`)

Interactive Jupyter notebooks with code, visualizations, and explanations.

#### `01_understanding_pyo3_bindings.ipynb`
**Topics**: PyO3, zero-copy, SIMD, memory layout, performance
**Duration**: ~45 minutes
**Level**: Intermediate

**You'll learn**:
- How PyO3 connects Python and Rust
- Zero-copy data sharing with NumPy
- SIMD vectorization benefits
- Memory layout and alignment
- Performance comparisons

**Prerequisites**: Basic Python, NumPy knowledge

---

#### `02_molecular_systems_and_energy.ipynb`
**Topics**: MD systems, states, energies, boxes, topology
**Duration**: ~60 minutes
**Level**: Intermediate

**You'll learn**:
- Creating molecular systems
- Initializing positions and velocities
- Energy tracking and components
- Periodic boundary conditions
- Complete system setup

**Prerequisites**: Basic molecular dynamics concepts

---

#### `03_performance_deep_dive.ipynb`
**Topics**: Benchmarking, profiling, optimization
**Duration**: ~90 minutes
**Level**: Advanced

**You'll learn**:
- Proper benchmarking techniques
- Memory profiling and layout
- SIMD performance analysis
- Cache effects and optimization
- Scaling behavior

**Prerequisites**: Performance analysis experience

### Scripts (`scripts/`)

Standalone Python scripts for learning and demonstration.

#### `explain_architecture.py`
**Type**: Educational walkthrough
**Duration**: ~20 minutes
**Level**: Beginner to Intermediate

**Covers**:
- Overall architecture (Polars pattern)
- PyO3 bridge explained
- Zero-copy data sharing
- SIMD acceleration
- Memory layout diagrams
- Performance benefits
- Typical workflow

**Usage**:
```bash
python scripts/explain_architecture.py
# Press Enter to navigate through sections
```

**Output**: Text-based interactive tutorial

---

#### `hands_on_demo.py`
**Type**: Interactive demonstrations
**Duration**: ~15 minutes
**Level**: Beginner

**Demos**:
1. Basic vector operations
2. NumPy integration
3. Performance comparison
4. Energy tracking
5. Simulation boxes
6. Complete system setup

**Usage**:
```bash
python scripts/hands_on_demo.py
# Follow the interactive prompts
```

**Output**: Live code execution with results

### Examples (`examples/`)

Production-quality example scripts.

| File | Topic | Level | Duration |
|------|-------|-------|----------|
| `01_basic_vectors.py` | Vector math + NumPy | Beginner | 5 min |
| `02_system_setup.py` | System initialization | Beginner | 10 min |
| `03_integrators.py` | MD integrators | Intermediate | 10 min |
| `04_advanced_sampling.py` | GaMD, EDS, REMD | Advanced | 15 min |
| `05_complete_workflow.py` | Full MD workflow | Intermediate | 20 min |

**Usage**:
```bash
cd examples
python 01_basic_vectors.py
python 02_system_setup.py
# etc.
```

## 🔧 Prerequisites

### Software Requirements
```bash
# Required
pip install numpy>=1.20
pip install maturin>=1.0

# Optional (for notebooks)
pip install jupyter matplotlib scipy pandas

# Optional (for development)
pip install pytest black ruff
```

### Build the Package
```bash
cd py-gromos
maturin develop --release
```

### Knowledge Requirements

| Material | Python | NumPy | Rust | MD | Performance |
|----------|--------|-------|------|-----|-------------|
| Scripts | ✓ | Basics | - | - | - |
| Notebook 01 | ✓ | ✓ | Basics | - | Basics |
| Notebook 02 | ✓ | ✓ | - | Basics | - |
| Notebook 03 | ✓ | ✓ | - | - | ✓ |
| Examples | ✓ | Basics | - | Basics | - |

## 🎓 Concepts Explained

### Core Concepts (Notebook 01)

**PyO3**: Rust ↔ Python bridge
```
Python code → PyO3 wrapper → Rust core
             (zero overhead)
```

**Zero-Copy**: Share memory, don't copy
```
NumPy array ──→ Same memory ←── Rust Vec
            (instant access)
```

**SIMD**: Process 4 values at once
```
Scalar: x + y requires 1 instruction per component
SIMD:   x + y requires 1 instruction total (4× faster!)
```

### MD Concepts (Notebook 02)

**State**: System snapshot
- Positions (x, y, z for each atom)
- Velocities (vx, vy, vz)
- Forces (fx, fy, fz)

**Energy**: Thermodynamic quantities
- Kinetic: ½mv²
- Potential: bonds + angles + LJ + Coulomb
- Total: conserved in NVE ensemble

**Box**: Periodic boundaries
- Rectangular: orthogonal box
- Triclinic: non-orthogonal (crystal systems)
- Vacuum: no periodicity

### Performance (Notebook 03)

**Benchmarking**: Measure correctly
- Warmup iterations
- Multiple runs
- Statistical analysis

**Memory**: Layout matters
- Contiguous: cache-friendly
- Aligned: SIMD-friendly
- Packed: space-efficient

**Optimization**:
- Batch operations
- Minimize copies
- Use float32
- Preallocate arrays

## 📊 Performance Summary

### Expected Speedups

| Operation | Pure Python | NumPy | GROMOS-RS | Best |
|-----------|------------|-------|-----------|------|
| Dot product | 1× | 20× | 66× | ✓ |
| Vector length | 1× | 15× | 72× | ✓ |
| Distance calc | 1× | 18× | 68× | ✓ |
| Array ops | 1× | 50× | 50× | ≈ |

### Why So Fast?

1. **Compiled Rust**: No interpreter overhead
2. **SIMD**: 4 operations at once
3. **Zero-copy**: No data movement
4. **Cache-friendly**: Contiguous memory
5. **Parallel**: Rayon threading (when applicable)

## 🚀 Next Steps

After completing the learning materials:

### For Users
1. Try the examples with your own data
2. Read the API reference for complete documentation
3. Check out the test suite for more examples
4. Join the community discussions

### For Developers
1. Read the Rust source code (`src/lib.rs`)
2. Study the PyO3 patterns used
3. Look at the gromos-rs core library
4. Consider contributing improvements

## 📝 Cheat Sheet

### Quick Reference

```python
import gromos
import numpy as np

# Vector operations (SIMD)
v1 = gromos.Vec3(1.0, 2.0, 3.0)
v2 = gromos.Vec3(4.0, 5.0, 6.0)
dist = v1.distance(v2)

# NumPy integration (zero-copy)
arr = v1.to_numpy()
v3 = gromos.Vec3.from_numpy(np.array([1, 2, 3], dtype=np.float32))

# System setup
state = gromos.State(num_atoms=1000, num_temp_groups=1, num_energy_groups=1)
box = gromos.Box.rectangular(5.0, 5.0, 5.0)
energy = gromos.Energy(num_temperature_groups=1, num_energy_groups=1)

# Data access
positions = np.random.rand(1000, 3).astype(np.float32)
state.set_positions(positions)
retrieved = state.positions()  # Zero-copy view

# Performance tip
for _ in range(1000):
    # Good: reuse arrays
    retrieved[:] = compute_new_positions()

    # Bad: allocate each time
    # retrieved = compute_new_positions()
```

## 🔗 References

### External Resources
- [Polars Project](https://github.com/pola-rs/polars) - Inspiration
- [PyO3 Documentation](https://pyo3.rs) - Python ↔ Rust bindings
- [NumPy Documentation](https://numpy.org/doc/) - Array operations
- [glam Documentation](https://docs.rs/glam/) - SIMD vectors

### Internal Documentation
- `README.md` - Installation and quick start
- `API_REFERENCE.md` - Complete API documentation
- Source code comments in `src/lib.rs`

## ❓ Troubleshooting

### Common Issues

**Q**: Import error - `ModuleNotFoundError: No module named 'gromos'`
**A**: Build the package first:
```bash
cd py-gromos
maturin develop --release
```

**Q**: Jupyter notebooks don't work
**A**: Install Jupyter:
```bash
pip install jupyter matplotlib
```

**Q**: Performance not as expected
**A**: Make sure you built with `--release` flag

**Q**: NumPy dtype mismatch
**A**: Use `np.float32` (not float64):
```python
arr = np.array([1, 2, 3], dtype=np.float32)  # ✓ Correct
arr = np.array([1, 2, 3])  # ✗ Wrong (float64)
```

## 💬 Feedback

Found an issue or have a suggestion?
- Open an issue on GitHub
- Contribute improvements via pull request
- Ask questions in discussions

## 📜 License

GPL-2.0 - Same as GROMOS

---

Happy learning! 🎉

For questions or issues, please refer to the main repository documentation.
