# Advanced Sampling Architecture for gromos-rs

## Overview

This document describes how REMD, EDS, and GaMD are integrated into gromos-rs following the GROMOS workflow architecture.

## GROMOS Workflow Pattern

```
mk_script (script generator)
    ↓
Multiple md/repex_mpi runs (actual simulation)
    ↓
GROMOS++ analysis tools (post-processing)
```

## Architecture Principles

1. **Separation of Concerns:**
   - `mk_script`: Generate job scripts, manage workflows
   - `md`/`repex_mpi`: Run simulations, write data files
   - Analysis tools: Process output files, compute statistics

2. **File-Based Communication:**
   - Input: `.imd` files contain simulation parameters
   - Output: `.cnf` (coordinates), `.trc` (trajectories), `.tre` (energies)
   - Special: `replica.dat` (REMD exchanges), `gamd.dat` (GaMD parameters)

3. **Restartability:**
   - Each run writes final coordinates for next run
   - Statistics files enable continuation
   - Scripts automatically chain together

## Implementation Strategy

### 1. REMD (Replica Exchange Molecular Dynamics)

#### MPI-Based Parallel Execution

**Binary:** `repex_mpi` (separate from main `md`)

**Usage:**
```bash
mpirun -np 8 repex_mpi \
    @topo system.top \
    @conf initial.cnf \
    @input remd.imd \
    @repdat replica.dat \
    @repout replica_%ID%.out
```

**Input File Block (remd.imd):**
```
REPLICA
# NRET  LRESCALE  RET(1:NRET)
  4     1         300.0 310.0 320.0 330.0
# TRIALS  EQUILIBRATE
  1000    100
END
```

**Output Files:**
- `replica.dat` - Exchange log (all replicas write here)
- `system_%ID%.cnf` - Per-replica coordinates
- `system_%ID%.trc` - Per-replica trajectories
- `system_%ID%.tre` - Per-replica energies
- `replica_%ID%.out` - Per-replica output logs

**Analysis Tool:** `rep_ana`
```bash
rep_ana @repdata replica.dat @input remd.imd
# Outputs: acceptance rates, temperature trajectories, diffusion
```

#### mk_script Integration

```bash
mk_script @sys peptide @bin repex_mpi @version mpi \
    @dir /work/remd @files
        topo     system.top
        input    remd.imd
        coord    equilibrated.cnf
        repdat   replica.dat
    @mpi 8  # Number of replicas
    @script 1 10  # 10 sequential runs
```

Generates: `peptide_1.run`, `peptide_2.run`, ..., `peptide_10.run`

Each script:
1. Runs MPI job with 8 processes
2. Reads previous `replica.dat` to continue exchanges
3. Writes updated `replica.dat` and per-replica files
4. Automatically submits next script

### 2. GaMD (Gaussian Accelerated MD)

#### Integrated into Main MD Binary

**Binary:** `md` (with GAMD block support)

**Three-Phase Workflow:**

**Phase 1: CMD Search (Classical MD statistics)**
```bash
md @topo system.top @conf initial.cnf @input gamd_cmd.imd @fin cmd_final.cnf
```

Input file `gamd_cmd.imd`:
```
GAMD
# SEARCH  FORM  THRESH
  1       2     1          # CMD search, total boost, lower bound
# SIGMA0_DIH  SIGMA0_TOT
  6.0         6.0
# EQUIL  WINDOW
  1000   0
END
```

**Phase 2: GaMD Search (Adaptive parameters)**
```bash
md @topo system.top @conf cmd_final.cnf @input gamd_search.imd @fin search_final.cnf
```

Input file `gamd_search.imd`:
```
GAMD
# SEARCH  FORM  THRESH
  2       2     1          # GaMD search, total boost, lower bound
# SIGMA0_DIH  SIGMA0_TOT
  6.0         6.0
# EQUIL  WINDOW
  0      1000              # 1000-step sliding window
END
```

**Phase 3: Production (Fixed parameters)**
```bash
md @topo system.top @conf search_final.cnf @input gamd_prod.imd @fin final.cnf
```

Input file `gamd_prod.imd`:
```
GAMD
# SEARCH  FORM  THRESH
  0       2     1          # Production, total boost, lower bound
# K_DIH    K_TOT      E_DIH    E_TOT
  0.0      0.001234   0.0      -1234.56
# EQUIL  WINDOW
  0      0
END
```

**Output Files:**
- Standard `.cnf`, `.trc`, `.tre` files
- `gamd_stats.dat` - Statistics from CMD/GaMD search
- `gamd_boost.dat` - Boost potential trajectory (for reweighting)

**Analysis Tool:** `gamd_ana`
```bash
gamd_ana @boost gamd_boost.dat @input gamd.imd
# Outputs: reweighted free energies, boost statistics
```

#### mk_script Integration

```bash
mk_script @sys peptide @bin md @dir /work/gamd \
    @joblist gamd_workflow.jobs @template mk_script.lib
```

Job file `gamd_workflow.jobs`:
```
JOBSCRIPTS
job_id  INPUT_FILE         RUN_AFTER  SUBDIR
1       gamd_cmd.imd       0          cmd
2       gamd_search.imd    1          search
3       gamd_prod.imd      2          prod
END
```

### 3. EDS (Enveloping Distribution Sampling)

#### Integrated into Main MD Binary

**Binary:** `md` (with EDS block support)

**Usage:**
```bash
md @topo system.top @conf initial.cnf @input eds.imd @fin final.cnf
```

**Input File Block (eds.imd):**
```
EDS
# NUMSTATES  FORM  S
  3          1     0.5          # 3 states, single-s, s=0.5
# E_OFFSETS (E_i^R for each state)
  0.0  10.0  20.0
# TEMP
  300.0
# SEARCH  E_MAX  E_MIN
  1       10.0   -50.0         # Enable AEDS parameter search
END
```

**Output Files:**
- Standard `.cnf`, `.trc`, `.tre` files
- `eds_stats.dat` - State visit statistics, round trips
- `eds_vr.dat` - Reference potential trajectory

**Analysis Tool:** `eds_ana`
```bash
eds_ana @vr eds_vr.dat @input eds.imd
# Outputs: state probabilities, free energies
```

## Implementation Tasks

### Core Modules (Keep as-is)
- [x] `src/gamd.rs` - GaMD logic
- [x] `src/eds.rs` - EDS logic
- [x] `src/remd.rs` - REMD logic
- [x] `src/replica.rs` - Replica management

### Input File Support (New)
- [ ] `src/io/input/gamd_block.rs` - Parse GAMD block
- [ ] `src/io/input/eds_block.rs` - Parse EDS block
- [ ] `src/io/input/replica_block.rs` - Parse REPLICA block

### Output File Writers (New)
- [ ] `src/io/output/repdat.rs` - Write replica.dat
- [ ] `src/io/output/gamd_stats.rs` - Write GaMD statistics
- [ ] `src/io/output/eds_stats.rs` - Write EDS statistics

### Binaries (Refactor)
- [ ] `src/bin/md.rs` - Add GAMD/EDS support via input blocks
- [ ] `src/bin/repex_mpi.rs` - MPI-based REMD (NEW)
- [ ] Remove `src/bin/remd.rs` (standalone - incorrect architecture)
- [ ] Remove `src/bin/gamd.rs` (standalone - incorrect architecture)
- [ ] Remove `src/bin/eds.rs` (standalone - incorrect architecture)

### Analysis Tools (New)
- [ ] `src/bin/rep_ana.rs` - Analyze replica.dat
- [ ] `src/bin/gamd_ana.rs` - Analyze GaMD boost potential
- [ ] `src/bin/eds_ana.rs` - Analyze EDS reference potential

### mk_script Enhancements
- [ ] `src/bin/mk_script.rs` - Add REMD support
- [ ] Add MPI job template support
- [ ] Add GaMD/EDS workflow templates
- [ ] Implement `@joblist` for sequential workflows

## File Format Specifications

### replica.dat Format
```
# REMD Exchange Data
# Replica_ID  Partner_ID  Run  λi    Ti      Epot_i      λj    Tj      Epot_j      Prob  Switch
1             2           1    0.0   300.0   -15234.5    0.0   310.0   -15198.3    0.876  1
2             1           1    0.0   310.0   -15198.3    0.0   300.0   -15234.5    0.876  1
3             4           1    0.0   320.0   -15165.2    0.0   330.0   -15142.8    0.654  0
4             3           1    0.0   330.0   -15142.8    0.0   320.0   -15165.2    0.654  0
```

### gamd_stats.dat Format
```
# GaMD Statistics
# Step  V_dih       V_tot       V_max_dih  V_min_dih  V_avg_dih  σ_dih   V_max_tot  V_min_tot  V_avg_tot  σ_tot    k_dih    k_tot      E_dih    E_tot
1000    -125.3      -15234.5    -100.2     -145.8     -123.4     12.3    -15100.0   -15350.0   -15220.0   58.4     0.0      0.001234   0.0      -15100.0
```

### eds_stats.dat Format
```
# EDS Statistics
# Step  State  V_i        V_R        E_offset   Visited  Round_trips
1000    1      -15234.5   -15220.1   0.0        1        0
1000    2      -15198.3   -15220.1   10.0       1        0
1000    3      -15165.2   -15220.1   20.0       0        0
```

## Migration Path

### Phase 1: Input File Support (Week 1)
1. Implement GAMD/EDS/REPLICA block parsers
2. Integrate into main `md` binary
3. Test with simple cases

### Phase 2: MPI REMD (Week 2)
1. Create `repex_mpi` binary
2. Implement replica.dat reading/writing
3. Test parallel execution

### Phase 3: Analysis Tools (Week 3)
1. Implement `rep_ana`
2. Implement `gamd_ana`
3. Implement `eds_ana`

### Phase 4: mk_script Integration (Week 4)
1. Add REMD/GaMD/EDS templates
2. Test complete workflows
3. Documentation and examples

## Benefits

1. **Matches GROMOS Architecture:** Users familiar with GROMOS can use gromos-rs the same way
2. **Scalability:** Long simulations broken into manageable chunks
3. **Restartability:** Can resume from any point
4. **Analysis Integration:** Works seamlessly with existing GROMOS++ tools
5. **Job Submission:** mk_script handles queue systems automatically

## Example Complete Workflow

```bash
# 1. Generate REMD workflow
mk_script @sys peptide @bin repex_mpi @version mpi \
    @dir /work/remd @mpi 8 @script 1 20 \
    @files topo system.top input remd.imd coord initial.cnf

# 2. Run first job (automatically chains to others)
./peptide_1.run

# 3. Wait for all 20 jobs to complete...

# 4. Analyze results
rep_ana @repdata replica.dat @input remd.imd

# 5. Extract specific replica trajectory
# (Using existing GROMOS++ tools)
```

This architecture enables gromos-rs to handle production-scale simulations while maintaining the flexibility and power of the GROMOS ecosystem.
