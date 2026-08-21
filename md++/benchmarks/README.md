# md++ performance benchmarks

Tracks md++ simulation performance across commits using
[airspeed velocity](https://asv.readthedocs.io) (asv).

The existing CI catches *correctness* regressions. This catches performance
ones -- a commit that makes the nonbonded loop 15% slower currently passes every
check we have.

## What is measured

Plain MD on three equilibrated protein-in-solvent systems, all using the same
reaction-field setup (0.8/1.4 nm cutoffs, SHAKE, pairlist every 5 steps), so
comparing across tiers isolates system size rather than method:

| Tier   | PDB    | Atoms   | ~ms/step (4 threads) |
|--------|--------|---------|----------------------|
| tiny   | `2ovn` | 11,904  | 8.2                  |
| medium | `1aki` | 45,530  | 33.9                 |
| large  | `1tua` | 172,541 | 136.8                |

md++ already reports its own wall times and a per-algorithm `TIMING` breakdown,
so each configuration runs *once* and many metrics are scraped from the output:

* `track_sim_walltime` -- simulation wall time (the headline number)
* `track_ms_per_step` -- wall time per step, comparable across systems
* `track_nonbonded`, `track_pairlist`, `track_nonbonded_shortrange`,
  `track_nonbonded_longrange`, `track_shake` -- per-algorithm breakdown
* `track_ns_per_day` -- throughput, informational only

The per-algorithm series are the point: they turn "this commit got slower" into
"the pairlist got slower".

Two deliberate choices worth knowing about:

* **Simulation time, not total time.** Initialisation is excluded. For the large
  system, reading the 29 MB input coordinates takes about four times longer than
  the simulation being measured.
* **ns/day is not the headline metric.** asv's regression detection assumes
  lower-is-better, so on a higher-is-better series every genuine speedup would be
  reported as a regression. Ignore asv's flags on that one series.

## Build variants

Selected by `GROMOS_VARIANT`; `scripts/asv_build.sh` maps each to its cmake
flags, so the combinations md++ forbids (OMP with MPI, CUKERNEL without OMP --
see `cmake/options.cmake`) cannot be expressed.

| Variant  | cmake flags               | Binary   | Launch                    |
|----------|---------------------------|----------|---------------------------|
| `serial` | *(none)*                  | `md`     | direct, 1 thread          |
| `omp`    | `-DOMP=ON`                | `md`     | direct, `OMP_NUM_THREADS` |
| `cuda`   | `-DOMP=ON -DCUKERNEL=ON`  | `md`     | direct + `INNERLOOP` CUDA |
| `mpi`    | `-DMPI=ON`                | `md_mpi` | `mpirun -np N`            |

Note that the CPU/GPU nonbonded kernel is a *runtime* choice made in the imd
`INNERLOOP` block, so the `cuda` variant is an OpenMP build that additionally
requests the CUDA kernel per run.

## Setup

Benchmark inputs are staged into a local cache once. They are deliberately not
committed: they total ~39 MB, and inputs must be byte-identical across every
commit in the timeline or the timeline stops measuring the code. The manifest
records a sha256 per file and is verified before every run.

```sh
python3 -m pip install --target /local/gromos/bench_tools asv virtualenv
python3 tools/stage_systems.py           # stage into /local/gromos/bench_systems
python3 tools/stage_systems.py --verify  # re-check against the manifest
```

The upstream protein set lives on shared storage and is read *only* by this
staging step. Everything else reads the local cache.

## Running

Use the `./asv` wrapper in this directory:

```sh
cd md++/benchmarks
./asv machine --yes                      # once per machine
./asv run HEAD^!                         # benchmark one commit
./asv run NEW                            # anything not yet benchmarked
./asv run 284700de..master --steps 20    # a sampled timeline
./asv publish && ./asv preview           # browsable results at :8080
./asv show                               # which commits have results
./asv compare <commit1> <commit2>        # side by side, flags regressions
```

The wrapper exists because asv is installed with `pip --target` into a directory
of its own, so it is not on `PATH` and needs `PYTHONPATH` set. It also cds to
its own directory so that asv picks up the `asv.conf.json` next to it, which
means it works when called by absolute path from anywhere.

If you would rather call `asv` directly:

```sh
export PYTHONPATH=/local/gromos/bench_tools
export PATH=/local/gromos/bench_tools/bin:$PATH
cd md++/benchmarks     # or pass --config .../asv.conf.json
```

A note if you ever run with `--config`: asv chdirs to the directory containing
the config file (`asv/main.py`) and resolves `repo`, `benchmark_dir` and the
`{conf_dir}` used by the build hooks relative to *that*, not to where you
invoked it. A config placed outside this directory therefore needs absolute
paths for all three, or the build hooks will not be found.

### Cost

Be deliberate here -- the full matrix is expensive:

| Scope | Time per commit |
|---|---|
| Full matrix (3 systems x 5 thread counts x 4 variants x 3 repeats) | **~2.7 h** |
| One variant, all systems | ~50 min |
| `tiny` only, all variants | ~38 min |
| `tiny` only, one variant, 1 repeat | ~3 min |

Two things drive this. Builds are one: roughly 230 translation units, so keep
ccache warm. The other is the large system's initialisation -- reading its 29 MB
input costs 54 s *per run*, which every repeat pays again; at 15 runs per variant
that is ~13 min of pure I/O before any physics happens.

So do not benchmark every commit. Use `./asv run NEW` or `--steps N` for a sampled
timeline, and narrow the scope for routine checks:

```sh
# fast smoke check
GROMOS_BENCH_TIERS=tiny GROMOS_BENCH_THREADS=4 GROMOS_BENCH_REPEATS=1 ./asv run HEAD^!

# nightly: one variant, full scaling curve
GROMOS_VARIANT=omp ./asv run NEW
```

### Machine-specific paths

`asv.conf.json` points `env_dir`, `results_dir` and `html_dir` at
`/local/gromos/bench_asv/...`, and the systems cache defaults to
`/local/gromos/bench_systems`. These are deliberately outside the repository --
asv checks project commits out into `env_dir`, and nesting that inside the
repository being checked out invites confusion -- but they *are* specific to
this machine. On another machine, either adjust the config or export
`GROMOS_BENCH_SYSTEMS` and `GROMOS_BENCH_SCRATCH`.

Expect `env_dir` to reach a few GB: each cached build is a ~55 MB binary and
`build_cache_size` keeps 8 of them.

### Prerequisites

* `nvcc` for the `cuda` variant. It is frequently not on `PATH`; set
  `GROMOS_NVCC` to its location.
* For low-noise results, set the CPU governor to `performance`
  (`sudo cpupower frequency-set -g performance`). `tools/calibrate.py` warns
  when it is not.

### History horizon

The suite supports commits back to `284700de` (2023-10-21), the
`gromosXX/` -> `md++/` rename. Older commits fail the build step with an
explicit message rather than a confusing cmake error.

## Environment variables

| Variable | Purpose |
|---|---|
| `GROMOS_VARIANT` | `serial`/`omp`/`cuda`/`mpi` (set by the asv matrix) |
| `GROMOS_BENCH_SYSTEMS` | staged systems cache |
| `GROMOS_BENCH_SYSTEMS_SRC` | upstream set, used only by `stage_systems.py` |
| `GROMOS_BENCH_SCRATCH` | run working directory; must be local disk |
| `GROMOS_BENCH_REPEATS` | repeats per configuration (default 3, best is kept) |
| `GROMOS_BENCH_TIERS` | restrict systems, e.g. `tiny` |
| `GROMOS_BENCH_THREADS` | restrict thread counts, e.g. `"1 4"` |
| `GROMOS_BENCH_PIN` | set `0` to disable `taskset` pinning |
| `GROMOS_BENCH_TRACE` | append one line per simulation to this file |
| `GROMOS_MD_BINARY` | use a prebuilt binary (standalone runs, bypasses asv) |
| `GROMOS_NVCC` | path to `nvcc` for the `cuda` variant |

## Standalone use

`tools/calibrate.py` runs the harness without asv -- useful for checking a local
build or re-tuning `NSTLIM` after a hardware change:

```sh
GROMOS_MD_BINARY=../BUILD/program/md python tools/calibrate.py --all
```

## Layout

```
asv                     wrapper: sets PYTHONPATH and runs from this directory
asv.conf.json           asv configuration and the build-variant matrix
scripts/asv_build.sh    cmake configure/build/install into asv's build cache
scripts/asv_install.sh  install a cached build into an asv environment
benchmarks/_imd.py      block-aware GROMOS imd reader/patcher
benchmarks/_omd.py      omd timing parser
benchmarks/_runner.py   variant dispatch, pinning, repeats
benchmarks/_systems.py  system registry and NSTLIM table
benchmarks/md_bench.py  the PlainMD benchmarks
tools/stage_systems.py  stage inputs into the local cache
tools/calibrate.py      measure per-step cost, suggest NSTLIM
```

## Notes on measurement noise

MD benchmarking is easy to get wrong. The harness disables all trajectory
writing in the imd (a production imd would otherwise have the benchmark
measuring the filesystem), runs in a local scratch directory, binds threads to
cores, pins the process to a fixed CPU set, and keeps the *minimum* across
repeats -- interference can only make a run slower, so the fastest observation
is closest to the machine's true capability.
