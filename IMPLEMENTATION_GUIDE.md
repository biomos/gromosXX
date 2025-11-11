# GROMOS Rust Port - Implementation Guide

## Quick Start

This guide walks you through the complete process of integrating Rust performance kernels into GROMOS.

## Phase 1: Build and Test Rust Library (30 minutes)

### Step 1: Build the Rust Library

```bash
cd gromos-rs

# Build in release mode with all optimizations
cargo build --release

# Verify the library was created
ls -lh target/release/libgromos_rs.a
```

### Step 2: Run Unit Tests

```bash
# Run all tests
cargo test --release

# Run specific test
cargo test --release test_lj_interaction

# Run with output
cargo test --release -- --nocapture
```

Expected output:
```
running 5 tests
test interaction::nonbonded::tests::test_innerloop_simple ... ok
test interaction::nonbonded::tests::test_lj_interaction ... ok
test interaction::nonbonded::tests::test_periodic_boundary ... ok
test math::tests::test_rectangular_boundary ... ok
test math::tests::test_vacuum_boundary ... ok

test result: ok. 5 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

### Step 3: Generate C Header

```bash
# Header is auto-generated during build
cat gromos_rs.h
```

## Phase 2: Benchmarking (1 hour)

### Step 1: Run Microbenchmarks

First, create a benchmark file:

```bash
mkdir -p benches
cat > benches/innerloop.rs << 'EOF'
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use gromos_rs::*;
use gromos_rs::interaction::nonbonded::*;
use gromos_rs::math::*;

fn bench_innerloop(c: &mut Criterion) {
    let mut group = c.benchmark_group("innerloop");

    for n_atoms in [100, 1000, 10000] {
        let positions: Vec<Vec3> = (0..n_atoms)
            .map(|i| Vec3::new(
                (i as f32 % 10.0),
                ((i / 10) as f32 % 10.0),
                ((i / 100) as f32 % 10.0),
            ))
            .collect();

        let charges: Vec<f32> = vec![0.5; n_atoms];
        let iac: Vec<u32> = vec![0; n_atoms];

        // Generate pairlist (approx 50 neighbors per atom)
        let mut pairlist = Vec::new();
        for i in 0..n_atoms.min(1000) {
            for j in (i + 1)..(i + 50).min(n_atoms) {
                pairlist.push((i as u32, j as u32));
            }
        }

        let lj_params = vec![vec![LJParameters { c6: 0.001, c12: 0.0001 }]];
        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431,
            crf_cut3i: 0.182215,
        };

        let periodicity = Vacuum;
        let mut storage = ForceStorage::new(n_atoms);

        group.bench_with_input(
            BenchmarkId::from_parameter(n_atoms),
            &n_atoms,
            |b, _| {
                b.iter(|| {
                    storage.clear();
                    lj_crf_innerloop(
                        black_box(&positions),
                        black_box(&charges),
                        black_box(&iac),
                        black_box(&pairlist),
                        black_box(&lj_params),
                        black_box(&crf),
                        black_box(&periodicity),
                        black_box(&mut storage),
                    );
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_innerloop);
criterion_main!(benches);
EOF
```

Run benchmarks:

```bash
# Run all benchmarks
cargo bench

# View HTML report
firefox target/criterion/report/index.html
```

### Step 2: Profile with Perf (Linux only)

```bash
# Build with debug symbols
cargo build --release

# Profile
perf record --call-graph dwarf ./target/release/examples/benchmark
perf report

# Generate flamegraph
cargo install flamegraph
cargo flamegraph --bench innerloop
```

## Phase 3: Integration with GROMOS (2-3 hours)

### Step 1: Modify CMake Build System

Edit `md++/CMakeLists.txt`:

```cmake
# Add after existing options
option(USE_RUST_KERNELS "Use Rust performance kernels" OFF)

if(USE_RUST_KERNELS)
    message(STATUS "Enabling Rust performance kernels")

    # Ensure Rust is available
    find_program(CARGO cargo REQUIRED)
    find_program(RUSTC rustc REQUIRED)

    # Build Rust library
    set(RUST_LIB_DIR ${CMAKE_SOURCE_DIR}/gromos-rs/target/release)
    set(RUST_LIB ${RUST_LIB_DIR}/libgromos_rs.a)

    add_custom_command(
        OUTPUT ${RUST_LIB}
        COMMAND cargo build --release
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/gromos-rs
        COMMENT "Building Rust library"
    )

    add_custom_target(gromos_rs_build DEPENDS ${RUST_LIB})

    # Import Rust library
    add_library(gromos_rs STATIC IMPORTED GLOBAL)
    set_target_properties(gromos_rs PROPERTIES
        IMPORTED_LOCATION ${RUST_LIB}
    )
    add_dependencies(gromos_rs gromos_rs_build)

    # Add include directory for header
    include_directories(${CMAKE_SOURCE_DIR}/gromos-rs)

    # Link with md++ programs
    target_link_libraries(libmd++ PUBLIC gromos_rs)

    # Platform-specific libraries
    if(UNIX AND NOT APPLE)
        target_link_libraries(libmd++ PUBLIC pthread dl m gcc_s)
    elseif(APPLE)
        target_link_libraries(libmd++ PUBLIC pthread)
    endif()

    # Define preprocessor macro
    target_compile_definitions(libmd++ PUBLIC USE_RUST_KERNELS)
endif()
```

### Step 2: Create C++ Wrapper

Create `md++/src/interaction/nonbonded/rust_wrapper.h`:

```cpp
#ifndef GROMOS_RUST_WRAPPER_H
#define GROMOS_RUST_WRAPPER_H

#ifdef USE_RUST_KERNELS

#include "../../../gromos-rs/gromos_rs.h"
#include "../../topology/topology.h"
#include "../../configuration/configuration.h"
#include "../interaction_types.h"

namespace interaction {

class RustNonbondedWrapper {
public:
    template<typename t_nonbonded_spec>
    static void calculate_nonbonded(
        topology::Topology& topo,
        configuration::Configuration& conf,
        const std::vector<std::pair<unsigned int, unsigned int>>& pairlist,
        const Nonbonded_Parameter& params,
        Storage& storage
    );
};

} // namespace interaction

#endif // USE_RUST_KERNELS
#endif // GROMOS_RUST_WRAPPER_H
```

### Step 3: Implement Wrapper

Create `md++/src/interaction/nonbonded/rust_wrapper.cc`:

```cpp
#ifdef USE_RUST_KERNELS

#include "rust_wrapper.h"
#include <vector>

namespace interaction {

template<typename t_nonbonded_spec>
void RustNonbondedWrapper::calculate_nonbonded(
    topology::Topology& topo,
    configuration::Configuration& conf,
    const std::vector<std::pair<unsigned int, unsigned int>>& pairlist,
    const Nonbonded_Parameter& params,
    Storage& storage
) {
    const unsigned int n_atoms = topo.num_atoms();
    const unsigned int n_pairs = pairlist.size();
    const unsigned int n_types = topo.num_atom_types();

    // Convert positions to flat array
    std::vector<float> positions_flat;
    positions_flat.reserve(n_atoms * 3);
    for (unsigned int i = 0; i < n_atoms; ++i) {
        const auto& pos = conf.current().pos(i);
        positions_flat.push_back(pos(0));
        positions_flat.push_back(pos(1));
        positions_flat.push_back(pos(2));
    }

    // Convert charges
    std::vector<float> charges_flat;
    charges_flat.reserve(n_atoms);
    for (unsigned int i = 0; i < n_atoms; ++i) {
        charges_flat.push_back(topo.charge()(i));
    }

    // Convert integer atom codes
    std::vector<uint32_t> iac_flat;
    iac_flat.reserve(n_atoms);
    for (unsigned int i = 0; i < n_atoms; ++i) {
        iac_flat.push_back(topo.iac(i));
    }

    // Convert pairlist
    std::vector<uint32_t> pairlist_flat;
    pairlist_flat.reserve(n_pairs * 2);
    for (const auto& pair : pairlist) {
        pairlist_flat.push_back(pair.first);
        pairlist_flat.push_back(pair.second);
    }

    // Convert LJ parameters
    std::vector<double> lj_params_flat;
    lj_params_flat.reserve(n_types * n_types * 2);
    for (unsigned int i = 0; i < n_types; ++i) {
        for (unsigned int j = 0; j < n_types; ++j) {
            const auto& lj = params.lj_parameter(i, j);
            lj_params_flat.push_back(lj.c6);
            lj_params_flat.push_back(lj.c12);
        }
    }

    // Get box information
    std::vector<double> box_data;
    uint32_t boundary_type;

    if constexpr (t_nonbonded_spec::boundary_type == math::vacuum) {
        boundary_type = 0;
        box_data = {0.0, 0.0, 0.0};
    } else if constexpr (t_nonbonded_spec::boundary_type == math::rectangular) {
        boundary_type = 1;
        const auto& box = conf.current().box;
        box_data = {box(0)(0), box(1)(1), box(2)(2)};
    } else {
        // Triclinic
        boundary_type = 2;
        // Not yet implemented
        throw std::runtime_error("Triclinic boundaries not yet supported in Rust kernel");
    }

    // CRF parameters
    double crf_cut = params.rf_cutoff;
    double crf_cut3 = crf_cut * crf_cut * crf_cut;
    double crf_2cut3i = 2.0 / crf_cut3;
    double crf_cut3i = 1.0 / crf_cut3;

    // Prepare output arrays
    std::vector<float> forces_flat(n_atoms * 3, 0.0f);
    std::vector<double> energies(2, 0.0);
    std::vector<double> virial(9, 0.0);

    // Call Rust kernel
    rust_lj_crf_innerloop(
        positions_flat.data(),
        charges_flat.data(),
        iac_flat.data(),
        n_atoms,
        pairlist_flat.data(),
        n_pairs,
        lj_params_flat.data(),
        n_types,
        box_data.data(),
        boundary_type,
        crf_cut,
        crf_2cut3i,
        crf_cut3i,
        forces_flat.data(),
        energies.data(),
        virial.data()
    );

    // Copy results back
    for (unsigned int i = 0; i < n_atoms; ++i) {
        storage.force(i)(0) += forces_flat[i * 3 + 0];
        storage.force(i)(1) += forces_flat[i * 3 + 1];
        storage.force(i)(2) += forces_flat[i * 3 + 2];
    }

    storage.energies.lj_energy += energies[0];
    storage.energies.crf_energy += energies[1];

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            storage.virial_tensor(i, j) += virial[i * 3 + j];
        }
    }
}

// Explicit template instantiations
// Add all used specializations here
template void RustNonbondedWrapper::calculate_nonbonded<...>(...);

} // namespace interaction

#endif // USE_RUST_KERNELS
```

### Step 4: Modify Nonbonded Calculation

Edit `md++/src/interaction/nonbonded/interaction/nonbonded_outerloop.cc`:

```cpp
#ifdef USE_RUST_KERNELS
#include "../rust_wrapper.h"
#endif

// In the main calculation function:
template<typename t_nonbonded_spec>
void Nonbonded_Outerloop<t_nonbonded_spec>::calculate(...) {

#ifdef USE_RUST_KERNELS
    // Use Rust kernel
    RustNonbondedWrapper::calculate_nonbonded<t_nonbonded_spec>(
        topo, conf, pairlist, m_param, storage
    );
#else
    // Use original C++ implementation
    for (const auto& pair : pairlist) {
        m_innerloop.lj_crf_innerloop(topo, conf, pair.first, pair.second,
                                     storage, periodicity);
    }
#endif
}
```

## Phase 4: Build and Test (1 hour)

### Step 1: Build GROMOS with Rust

```bash
cd md++
mkdir build-rust
cd build-rust

# Configure with Rust enabled
cmake .. -DUSE_RUST_KERNELS=ON -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)

# Verify Rust library was linked
ldd program/md | grep gromos_rs
```

### Step 2: Run Test Simulations

```bash
# Small test (argon gas)
./program/md << EOF
@topo argon.top
@conf argon.cnf
@input md.inp
@fin final.cnf
@trc trajectory.trc
EOF

# Compare energy conservation with C++ version
tail -n 100 md.out | grep "total energy"
```

### Step 3: Validation Tests

```bash
# Energy conservation test (NVE ensemble)
# Should have drift < 1e-4 kJ/mol per ns

# Force comparison test
# Run same input with C++ and Rust, compare forces

# Long simulation test
# Run 1 ns simulation, check RMSD from crystal structure
```

## Phase 5: Performance Validation (1 hour)

### Benchmark Suite

Create `benchmark_suite.sh`:

```bash
#!/bin/bash

SYSTEMS=("lysozyme" "dhfr" "membrane")
STEPS=10000

for system in "${SYSTEMS[@]}"; do
    echo "Benchmarking $system"

    # C++ version
    echo "  C++ version..."
    time ./md_cpp < ${system}.inp > ${system}_cpp.out

    # Rust version
    echo "  Rust version..."
    time ./md_rust < ${system}.inp > ${system}_rust.out

    # Extract timings
    cpp_time=$(grep "Total time" ${system}_cpp.out | awk '{print $3}')
    rust_time=$(grep "Total time" ${system}_rust.out | awk '{print $3}')

    speedup=$(echo "scale=2; $cpp_time / $rust_time" | bc)
    echo "  Speedup: ${speedup}x"
done
```

Expected results:
- Lysozyme (2,600 atoms): 1.8-2.2x speedup
- DHFR (23,500 atoms): 2.0-2.5x speedup
- Membrane (85,000 atoms): 2.3-2.8x speedup

## Troubleshooting

### Common Issues

#### 1. Linking Errors

```
undefined reference to 'rust_lj_crf_innerloop'
```

**Solution**: Ensure library order is correct:
```cmake
target_link_libraries(libmd++ PUBLIC gromos_rs pthread dl m)
```

#### 2. Different Results

```
Energy drift different from C++ version
```

**Solution**: Check floating-point precision:
- Ensure same CRF parameters
- Verify pairlist generation identical
- Check periodic boundary wrapping

#### 3. Performance Not Improving

```
Rust version slower than C++
```

**Solution**: Verify optimization flags:
```bash
# Rebuild with native CPU features
cd gromos-rs
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

#### 4. Segmentation Fault

```
Segfault in rust_lj_crf_innerloop
```

**Solution**: Check array sizes:
- Verify `n_atoms`, `n_pairs`, `n_types` are correct
- Ensure all pointers are non-null
- Check that arrays are large enough

## Performance Tuning

### CPU-Specific Optimizations

```bash
# Intel Skylake+
RUSTFLAGS="-C target-cpu=skylake" cargo build --release

# AMD Zen 2+
RUSTFLAGS="-C target-cpu=znver2" cargo build --release

# Apple Silicon
RUSTFLAGS="-C target-cpu=apple-m1" cargo build --release
```

### Thread Tuning

```bash
# Set thread count
export RAYON_NUM_THREADS=16

# Disable parallel for small systems
export RAYON_NUM_THREADS=1  # For < 1000 atoms
```

### Profile-Guided Optimization (Advanced)

```bash
# 1. Instrument build
cd gromos-rs
RUSTFLAGS="-Cprofile-generate=/tmp/pgo-data" cargo build --release

# 2. Run representative workload
cd ../md++/build-rust
./program/md < typical_system.inp

# 3. Optimized build
cd ../../gromos-rs
RUSTFLAGS="-Cprofile-use=/tmp/pgo-data" cargo build --release

# Expected: Additional 5-10% speedup
```

## Success Criteria

✅ All unit tests pass
✅ Energy conservation test passes (drift < 1e-4 kJ/mol/ns)
✅ Forces match C++ version to within 1e-6
✅ 2x speedup on representative benchmark
✅ No segfaults or memory leaks in 1 µs simulation
✅ Bit-reproducible results across runs

## Next Steps

After successful integration:

1. **Port pairlist generation** to Rust
2. **Add SIMD vectorization** for distance calculations
3. **Implement triclinic boundaries**
4. **Port SHAKE/SETTLE** constraints
5. **Add GPU offload** support

## Support

- GitHub Issues: https://github.com/gromos/gromos-rs/issues
- Email: biomos@igc.phys.chem.ethz.ch
- Documentation: See RUST_PORT_STRATEGY.md
