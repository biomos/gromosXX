# GROMOS Rust Validation Plan

## Problem Statement

We have translated GROMOS algorithms to Rust, but **haven't validated** that:
1. It can read GROMOS input files (topology, coordinates, parameters)
2. It produces the same output as C++ GROMOS
3. Results are bit-exact or numerically equivalent

## Validation Strategy

### Phase 1: File I/O (Read GROMOS Formats)

GROMOS uses specific file formats that we need to parse:

#### Input Files
- **`.top`** - Topology file (atoms, bonds, force field parameters)
- **`.cnf`** - Configuration file (positions, velocities, box)
- **`.imd`** - Input parameters (simulation settings)
- **`.ptp`** - Parameter topology
- **`.ifp`** - Force field parameters

#### Output Files
- **`.cnf`** - Final configuration
- **`.trc`** - Trajectory (positions over time)
- **`.tre`** - Energy trajectory
- **`.trv`** - Velocity trajectory

### Phase 2: Validation Tests

#### Level 1: Single-Point Energy
```
Input:  positions.cnf, topology.top
Action: Calculate forces and energies (no MD)
Check:  Rust energies == C++ energies (within tolerance)
```

#### Level 2: Short MD Trajectory
```
Input:  start.cnf, topology.top, md.imd
Action: Run 100 steps of MD
Check:  Rust trajectory ≈ C++ trajectory
```

#### Level 3: Long Simulation
```
Input:  system.cnf, topology.top, production.imd
Action: Run 10,000 steps
Check:  Energy conservation, temperature, pressure
```

### Phase 3: Automated Comparison

Create test suite that:
1. Runs same input through C++ and Rust
2. Compares outputs automatically
3. Reports any differences

## Implementation Plan

### Step 1: Create I/O Module

```rust
// gromos-rs/src/io/mod.rs
pub mod topology;  // Read .top files
pub mod coordinate; // Read .cnf files
pub mod trajectory; // Write .trc files
pub mod energy;    // Write .tre files
```

### Step 2: Validation Test Framework

```rust
// gromos-rs/tests/validation.rs
#[test]
fn test_single_point_energy() {
    // Load reference data
    let ref_energy = load_cpp_reference("test_data/ref_energy.dat");

    // Calculate with Rust
    let topo = read_topology("test_data/system.top");
    let conf = read_coordinate("test_data/system.cnf");
    let rust_energy = calculate_energy(&topo, &conf);

    // Compare
    assert_relative_eq!(rust_energy, ref_energy, epsilon = 1e-6);
}
```

### Step 3: Integration Tests

```bash
# Run C++ GROMOS
cd test_data/argon_gas
../../md++/program/md < argon.inp > cpp_output.txt

# Run Rust GROMOS
cargo run --release --example md_argon > rust_output.txt

# Compare results
python compare_outputs.py cpp_output.txt rust_output.txt
```

## What's Needed Now

### Immediate Priority: I/O Modules

1. **Topology Parser** - Read GROMOS `.top` files
2. **Coordinate Parser** - Read GROMOS `.cnf` files
3. **Trajectory Writer** - Write `.trc` files compatible with GROMOS tools

### Test Data Requirements

Need from C++ GROMOS:
- Small test systems (argon gas, water box, small protein)
- Reference energies and forces
- Reference trajectories

### Validation Tools

1. **Energy comparison script**
2. **Trajectory RMSD calculator**
3. **Statistical analysis** (mean, std dev of energies)

## Timeline

### Week 1: I/O Foundation
- [ ] Parse GROMOS topology format
- [ ] Parse GROMOS coordinate format
- [ ] Read force field parameters

### Week 2: Single-Point Validation
- [ ] Calculate energies for test systems
- [ ] Compare with C++ reference
- [ ] Debug any discrepancies

### Week 3: MD Validation
- [ ] Run short MD simulations
- [ ] Compare trajectories
- [ ] Validate energy conservation

### Week 4: Production Testing
- [ ] Large system tests
- [ ] Long time simulations
- [ ] Performance benchmarks

## Success Criteria

### Must Have
- ✅ Reads standard GROMOS input files
- ✅ Single-point energies within 1e-6 of C++
- ✅ 100-step trajectory within 1e-3 RMSD
- ✅ Energy conservation < 1e-4 kJ/mol/ns

### Should Have
- ✅ Automatic validation test suite
- ✅ Comparison tools for all outputs
- ✅ Documentation of differences

### Nice to Have
- ✅ Bit-exact results (if possible)
- ✅ Validation on diverse systems
- ✅ Continuous validation in CI/CD

## Known Challenges

### Numerical Precision
- Floating-point operations may differ slightly
- Acceptable tolerance: 1e-6 to 1e-4 depending on quantity

### Random Number Generation
- Stochastic methods (Langevin, etc.) need same RNG
- May need to implement GROMOS RNG exactly

### Compiler Differences
- LLVM (Rust) vs GCC/Intel (C++) may give slightly different results
- This is acceptable if within physical tolerance

## Next Steps

1. **Create I/O module structure** (today)
2. **Implement topology parser** (this week)
3. **Add test data** (small systems)
4. **Run first validation** (single-point energy)
5. **Iterate on differences**

---

**Current Status**: ⚠️ Translation complete but NOT validated
**Target Status**: ✅ Validated against C++ GROMOS reference
**ETA**: 2-4 weeks for full validation
