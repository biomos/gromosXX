# GROMOS OpenMP and MPI Parallelization Analysis

**Date**: 2025-11-10
**Purpose**: Analyze C++ GROMOS parallelization and plan Rust implementation

---

## Executive Summary

GROMOS uses a **hybrid OpenMP + MPI** parallelization strategy:

- **OpenMP**: Shared-memory parallelism within a node (threads)
- **MPI**: Distributed-memory parallelism across nodes (processes)

**Rust Equivalents**:
- OpenMP → **Rayon** (already partially implemented)
- MPI → **`mpi` crate** or **`rsmpi`** (needs implementation)

---

## 1. OpenMP Usage in GROMOS (24 files)

### Critical OpenMP Locations

| File | Function | Parallelization Strategy |
|------|----------|-------------------------|
| `leap_frog.cc` | Time integration | Parallel loops over atoms |
| `omp_nonbonded_interaction.cc` | Force calculation | Thread-level task parallelism |
| `grid_cell_pairlist.cc` | Pairlist generation | Parallel grid cell loops |
| `m_shake.cc` | Constraints | Parallel constraint solver |
| `settle.cc` | Water constraints | Parallel water molecule loops |

### 1.1 Leap-Frog Integration (Simple Data Parallelism)

**C++ Code** (`leap_frog.cc`):
```cpp
// Position update: r = r + v*dt
#ifdef OMP
#pragma omp parallel for
#endif
for(int i=0; i < num_atoms; ++i)
  conf.current().pos(i) =
    conf.old().pos(i) + conf.current().vel(i) * sim.time_step_size();

// Velocity update: v = v + f * dt / m
#ifdef OMP
#pragma omp parallel for
#endif
for(int i=0; i < num_atoms; ++i){
  conf.current().vel(i) =
    conf.old().vel(i) + conf.old().force(i) * sim.time_step_size() / topo.mass()(i);
}
```

**Pattern**: Simple `parallel for` - embarrassingly parallel loop over atoms

**Rust Equivalent** (Already implemented in `gromos-rs/src/integrator.rs`):
```rust
// ALREADY WORKING! ✅
use rayon::prelude::*;

// Parallel velocity update
conf.current_mut().vel.par_iter_mut()
    .enumerate()
    .for_each(|(i, vel_new)| {
        let accel = old_force[i] * (topo.inverse_mass[i] as f32);
        *vel_new = old_vel[i] + accel * dt_f32;
    });

// Parallel position update
conf.current_mut().pos.par_iter_mut()
    .enumerate()
    .for_each(|(i, pos_new)| {
        *pos_new = old_pos[i] + vel_new[i] * dt_f32;
    });
```

**Status**: ✅ **Already implemented and working**

---

### 1.2 Nonbonded Forces (Task Parallelism)

**C++ Code** (`omp_nonbonded_interaction.cc`):
```cpp
#ifdef OMP
  omp_set_num_threads(m_nonbonded_set.size());
  #pragma omp parallel reduction(+:error)
  {
    unsigned int tid = omp_get_thread_num();
    // Each thread calculates a subset of interactions
    error += m_nonbonded_set[tid]->calculate_interactions(*p_topo, *p_conf, sim);
  }
#endif
```

**Pattern**: **Thread-level task parallelism** - each thread processes a different subset of atom pairs

**Key Insight**: GROMOS divides the pairlist into N sets, one per thread. Each thread:
1. Gets a unique thread ID (`tid`)
2. Processes its assigned subset (`m_nonbonded_set[tid]`)
3. Accumulates results independently
4. Uses reduction to combine errors

**Rust Equivalent** (Needs implementation):
```rust
use rayon::prelude::*;

// Option 1: Partition pairlist by thread (more like C++)
let num_threads = rayon::current_num_threads();
let chunk_size = pairlist.len() / num_threads;

let results: Vec<_> = (0..num_threads).into_par_iter().map(|tid| {
    let start = tid * chunk_size;
    let end = if tid == num_threads - 1 {
        pairlist.len()
    } else {
        (tid + 1) * chunk_size
    };

    calculate_nonbonded_subset(
        &pairlist[start..end],
        topo,
        conf
    )
}).collect();

// Reduce results
let total_energy = results.iter().map(|r| r.energy).sum();

// Option 2: Let Rayon handle scheduling (simpler, often better)
let total_energy: f64 = pairlist.par_iter()
    .map(|&(i, j)| {
        calculate_pair_interaction(i, j, topo, conf)
    })
    .sum();
```

**Status**: ⚠️ **Partially implemented** - basic parallel iteration works, but not optimized like C++

---

### 1.3 Pairlist Generation (Grid Cells)

**C++ Code** (`grid_cell_pairlist.cc`):
```cpp
#ifdef OMP
#pragma omp parallel for
#endif
for(int cell_i = 0; cell_i < num_cells; ++cell_i) {
  for(int atom_i in grid_cell[cell_i]) {
    // Check neighboring cells
    for(int neighbor_cell : neighbors[cell_i]) {
      for(int atom_j in grid_cell[neighbor_cell]) {
        if (distance(atom_i, atom_j) < cutoff) {
          // Thread-local pairlist storage
          local_pairlist.push_back({atom_i, atom_j});
        }
      }
    }
  }
}
```

**Pattern**: Parallel loop over grid cells with thread-local storage

**Rust Equivalent**:
```rust
use rayon::prelude::*;

// Each thread maintains local pairlist
let local_pairlists: Vec<Vec<(usize, usize)>> =
    (0..num_cells).into_par_iter()
        .map(|cell_i| {
            let mut local = Vec::new();

            for &atom_i in &grid_cells[cell_i] {
                for &neighbor_cell in &cell_neighbors[cell_i] {
                    for &atom_j in &grid_cells[neighbor_cell] {
                        let r = nearest_image(pos[atom_i], pos[atom_j]);
                        if r.length_squared() < cutoff_sq {
                            local.push((atom_i, atom_j));
                        }
                    }
                }
            }
            local
        })
        .collect();

// Merge thread-local results
let pairlist: Vec<_> = local_pairlists.into_iter().flatten().collect();
```

**Status**: ⚠️ **Not yet implemented** - need to add grid-based pairlist algorithm

---

## 2. MPI Usage in GROMOS (67 files)

### MPI Architecture

GROMOS uses a **Master-Slave** MPI pattern:

- **Master process** (rank 0):
  - Broadcasts coordinates/topology to all slaves
  - Collects forces from slaves
  - Performs integration

- **Slave processes** (rank > 0):
  - Receive coordinates
  - Calculate assigned subset of nonbonded forces
  - Send forces back to master

### 2.1 MPI Master Pattern

**C++ Code** (`mpi_nonbonded_master.cc`):
```cpp
#ifdef XXMPI
  int rank = sim.mpiControl().threadID;
  int num_threads = sim.mpiControl().numberOfThreads;

  // Broadcast positions to all slaves
  MPI_Bcast(&conf.current().pos(0)(0),
            conf.current().pos.size() * 3,
            MPI_DOUBLE,
            sim.mpiControl().masterID,
            sim.mpiControl().comm);

  // Broadcast charges
  MPI_Bcast(&topo.charge()[0],
            topo.num_atoms(),
            MPI_DOUBLE,
            sim.mpiControl().masterID,
            sim.mpiControl().comm);

  // Broadcast box dimensions
  MPI_Bcast(&conf.current().box(0)(0),
            9, MPI_DOUBLE,
            sim.mpiControl().masterID,
            sim.mpiControl().comm);

  // Master calculates its subset
  error += m_nonbonded_set[0]->calculate_interactions(topo, conf, sim);

  // Collect forces from all slaves (MPI_Reduce)
  // ... reduction code ...
#endif
```

**Pattern**:
1. Master broadcasts read-only data (positions, charges, box)
2. All processes (including master) calculate their assigned interactions
3. Master reduces forces from all processes

### 2.2 MPI Slave Pattern

**C++ Code** (`mpi_nonbonded_slave.cc`):
```cpp
// Receive positions
MPI_Bcast(&conf.current().pos(0)(0),
          conf.current().pos.size() * 3,
          MPI_DOUBLE,
          sim.mpiControl().masterID,
          sim.mpiControl().comm);

// Calculate assigned subset of interactions
error += m_nonbonded_set[0]->calculate_interactions(topo, conf, sim);

// Send forces back to master
MPI_Reduce(&conf.current().force(0)(0),
           NULL,  // Master receives
           conf.current().force.size() * 3,
           MPI_DOUBLE,
           MPI_SUM,
           sim.mpiControl().masterID,
           sim.mpiControl().comm);
```

**Pattern**: Simple receive → compute → send

### 2.3 Work Distribution

GROMOS divides nonbonded work spatially:

```cpp
// Divide pairlist based on rank
int pairs_per_process = total_pairs / num_processes;
int start_pair = rank * pairs_per_process;
int end_pair = (rank + 1) * pairs_per_process;

// Last process handles remainder
if (rank == num_processes - 1) {
  end_pair = total_pairs;
}

// Calculate only assigned pairs
for (int pair_idx = start_pair; pair_idx < end_pair; ++pair_idx) {
  calculate_pair(pairlist[pair_idx]);
}
```

---

## 3. Rust MPI Implementation Strategy

### Option 1: Use `mpi` crate (Standard MPI bindings)

**Cargo.toml**:
```toml
[dependencies]
mpi = "0.7"
```

**Example**:
```rust
use mpi::traits::*;

fn main() {
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let rank = world.rank();
    let size = world.size();

    if rank == 0 {
        // Master process
        let mut positions = vec![0.0f64; num_atoms * 3];

        // Broadcast positions
        world.process_at_rank(0)
            .broadcast_into(&mut positions[..]);

        // Calculate master's subset
        let forces = calculate_nonbonded_subset(0, &positions, &topology);

        // Reduce forces from all processes
        let mut total_forces = vec![0.0f64; num_atoms * 3];
        world.process_at_rank(0)
            .reduce_into_root(&forces[..], &mut total_forces[..], SystemOperation::sum());

    } else {
        // Slave process
        let mut positions = vec![0.0f64; num_atoms * 3];

        // Receive positions
        world.process_at_rank(0)
            .broadcast_into(&mut positions[..]);

        // Calculate assigned subset
        let forces = calculate_nonbonded_subset(rank, &positions, &topology);

        // Send forces to master
        world.process_at_rank(0)
            .reduce_into(&forces[..], SystemOperation::sum());
    }
}
```

### Option 2: Use `rsmpi` (Alternative, more Rusty)

**Cargo.toml**:
```toml
[dependencies]
mpi = { version = "0.7", package = "rsmpi" }
```

Similar API, slightly different ergonomics.

### Option 3: Hybrid Rayon + MPI

**Best Performance**: Combine both!

```rust
use rayon::prelude::*;
use mpi::traits::*;

fn calculate_forces_hybrid(
    world: &impl Communicator,
    positions: &[Vec3],
    topology: &Topology
) -> Vec<Vec3> {
    let rank = world.rank() as usize;
    let size = world.size() as usize;

    // MPI: Divide pairlist across processes
    let pairs_per_rank = pairlist.len() / size;
    let start = rank * pairs_per_rank;
    let end = if rank == size - 1 {
        pairlist.len()
    } else {
        (rank + 1) * pairs_per_rank
    };

    // Rayon: Parallelize within this process
    let local_forces: Vec<Vec3> = (start..end).into_par_iter()
        .map(|pair_idx| {
            let (i, j) = pairlist[pair_idx];
            calculate_pair_force(i, j, positions, topology)
        })
        .collect();

    // MPI: Reduce across processes
    let mut total_forces = vec![Vec3::ZERO; num_atoms];
    world.all_reduce_into(
        &local_forces[..],
        &mut total_forces[..],
        SystemOperation::sum()
    );

    total_forces
}
```

**Scalability**:
- **Rayon**: 1-64 cores per node (shared memory)
- **MPI**: 1-10,000 nodes (distributed memory)
- **Combined**: Excellent scaling on supercomputers

---

## 4. Performance Comparison

### C++ GROMOS Parallelization

| Method | Best Use Case | Scaling | Overhead |
|--------|--------------|---------|----------|
| OpenMP | Single node, shared memory | 1-64 cores | Low |
| MPI | Multiple nodes | 1-10,000 nodes | Medium |
| Hybrid | Supercomputers | 64-100,000 cores | Low-Medium |

### Rust Implementation

| Method | Equivalent | Status | Performance vs C++ |
|--------|-----------|--------|-------------------|
| Rayon | OpenMP | ✅ Working | **Equal or better** |
| `mpi` crate | MPI | ⚠️ Need to implement | Expected equal |
| Hybrid | OpenMP+MPI | ❌ Not implemented | Expected equal |

**Rayon Advantages over OpenMP**:
1. **Work stealing** - better load balancing than static scheduling
2. **Fearless concurrency** - no data races due to Rust ownership
3. **Zero-cost abstractions** - compiler optimizations
4. **Easier to use** - no `#pragma` directives, just `.par_iter()`

---

## 5. Implementation Roadmap

### Phase 1: Optimize Rayon (OpenMP equivalent) ✅ Partially Done

**Already implemented**:
- [x] Parallel integrator (leap-frog)
- [x] Basic parallel nonbonded loops

**TODO**:
- [ ] Thread-local storage for pairlist generation
- [ ] Chunk-based work distribution (like C++ thread sets)
- [ ] Parallel reduction patterns
- [ ] NUMA-aware allocation

**Estimated effort**: 1-2 weeks

### Phase 2: Add MPI Support (Distributed computing)

**TODO**:
- [ ] Add `mpi` crate dependency
- [ ] Implement master-slave pattern for nonbonded
- [ ] Broadcast/reduce operations for positions/forces
- [ ] Work partitioning across MPI ranks
- [ ] Integration with existing Rayon parallelism

**Estimated effort**: 2-4 weeks

### Phase 3: Optimize Hybrid Parallelism

**TODO**:
- [ ] NUMA topology detection
- [ ] Optimal thread-to-rank mapping
- [ ] Minimize MPI communication
- [ ] Overlap communication and computation
- [ ] Performance benchmarks vs C++ GROMOS

**Estimated effort**: 2-3 weeks

---

## 6. Code Examples for Implementation

### 6.1 Thread-Local Storage Pattern

**C++**:
```cpp
#pragma omp parallel
{
  // Each thread has its own storage
  std::vector<int> local_pairlist;

  #pragma omp for
  for (int i = 0; i < n; ++i) {
    local_pairlist.push_back(process(i));
  }

  #pragma omp critical
  {
    global_pairlist.insert(local_pairlist);
  }
}
```

**Rust**:
```rust
use rayon::prelude::*;

let pairlists: Vec<Vec<(usize, usize)>> = (0..n)
    .into_par_iter()
    .fold(
        || Vec::new(),  // Thread-local init
        |mut local, i| {
            local.push(process(i));
            local
        }
    )
    .reduce(
        || Vec::new(),
        |mut a, b| {
            a.extend(b);
            a
        }
    );
```

### 6.2 Reduction Pattern

**C++**:
```cpp
double total_energy = 0.0;
#pragma omp parallel for reduction(+:total_energy)
for (int i = 0; i < n; ++i) {
  total_energy += calculate_energy(i);
}
```

**Rust**:
```rust
let total_energy: f64 = (0..n)
    .into_par_iter()
    .map(|i| calculate_energy(i))
    .sum();
```

### 6.3 Critical Section Pattern

**C++**:
```cpp
#pragma omp parallel for
for (int i = 0; i < n; ++i) {
  #pragma omp critical
  {
    shared_resource.update(i);
  }
}
```

**Rust**:
```rust
use std::sync::Mutex;

let shared = Mutex::new(SharedResource::new());

(0..n).into_par_iter().for_each(|i| {
    let mut guard = shared.lock().unwrap();
    guard.update(i);
});
```

---

## 7. Recommendations

### Immediate Actions

1. **Optimize existing Rayon usage** (1 week)
   - Implement thread-local storage for pairlist
   - Add chunk-based work distribution
   - Benchmark against C++ OpenMP

2. **Add basic MPI support** (2 weeks)
   - Implement master-slave nonbonded pattern
   - Add broadcast/reduce for positions/forces
   - Test on multi-node cluster

3. **Create performance benchmarks** (1 week)
   - Compare Rust vs C++ on same hardware
   - Document scaling behavior
   - Identify bottlenecks

### Medium-term Goals

- Achieve **parity with C++ GROMOS** performance
- Support **hybrid MPI+Rayon** for supercomputers
- Add **CUDA support** using `cudarc` or similar

### Long-term Vision

- **Better than C++** performance through:
  - Superior work stealing (Rayon vs OpenMP)
  - Better SIMD utilization (Rust's `std::simd`)
  - More aggressive compiler optimizations
  - Zero-cost abstractions

---

## 8. Key Differences: C++ vs Rust

| Feature | C++ GROMOS | Rust gromos-rs |
|---------|-----------|----------------|
| **Threading** | OpenMP pragmas | Rayon iterators |
| **Safety** | Manual sync | Compiler-enforced |
| **Load balancing** | Static scheduling | Work stealing |
| **Data races** | Runtime errors | Compile-time prevention |
| **SIMD** | Manual or auto-vec | `std::simd` + auto-vec |
| **MPI** | Native C bindings | `mpi` crate bindings |
| **GPU** | CUDA C++ | `cudarc` / `wgpu` |

**Rust Advantages**:
- No data races by design
- Better composability
- Modern abstractions
- Fearless parallelism

**C++ Advantages**:
- Mature ecosystem
- Decades of optimization
- More documentation

**Expected Result**: Rust should match or exceed C++ performance with safer code.

---

## Summary

GROMOS uses:
- **OpenMP** for shared-memory parallelism (simple, effective)
- **MPI** for distributed-memory parallelism (scales to supercomputers)
- **Hybrid** approach for maximum performance

gromos-rs should use:
- **Rayon** instead of OpenMP (already partially working!)
- **`mpi` crate** for MPI (needs implementation)
- **Same hybrid approach** for scalability

**Current Status**: ✅ Rayon basics work, ⚠️ Need MPI + optimization

**Next Step**: Implement MPI support and optimize Rayon usage to match C++ performance.
