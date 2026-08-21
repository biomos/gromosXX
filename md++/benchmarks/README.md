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
* `track_init_walltime` -- setup before step one; scales far worse than
  the simulation does, so worth watching separately
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

### Timelines

A timeline over the full matrix builds md++ four times per commit, which is
rarely what routine tracking needs. `quick.conf.json` is the same
configuration restricted to the `omp` variant -- asv's `-E` flag selects only
the environment type and Python version, so a separate config is the only way to
pick a single build variant.

```sh
# five commits sampled from the last eight on master, one variant
GROMOS_BENCH_TIERS=tiny GROMOS_BENCH_THREADS=4 \
  ./asv run --config quick.conf.json --bench "PlainMD.track_" \
            --steps 5 master~7..master

# the normal nightly invocation: whatever is not yet benchmarked
./asv run --config quick.conf.json --skip-existing-successful NEW

./asv publish && ./asv preview
```

Add `--skip-existing-successful` to repeat invocations so commits already
measured are not redone. Leave `GROMOS_BENCH_REPEATS` at its default of 3 for
real tracking: a single repeat leaves a run-to-run spread of roughly 1.4% on
this machine, which is the floor on what a regression has to exceed to be
visible.

### Cost

All figures below are measured on this machine, not estimated.

| Scope (per commit) | Time |
|---|---|
| Full matrix: 4 variants, 3 systems, 5 thread counts, 3 repeats | **2.8 h** |
| One variant, 3 systems, 5 thread counts, 3 repeats | 49 min |
| One variant, `tiny`+`medium`, 5 thread counts, 3 repeats | 25 min |
| One variant, `tiny`+`medium`, 5 thread counts, 2 repeats | 17 min |
| One variant, `tiny` only, 5 thread counts, 3 repeats | 12 min |

Where the time actually goes, for the full matrix:

| Stage | Time | Share |
|---|---|---|
| Builds (4 x 65 s, ccache warm) | 4.2 min | 2.5% |
| Initialisation | 46.9 min | 28.3% |
| Simulation -- the part being measured | 114.8 min | 69.2% |

**Builds are not the bottleneck.** They dominate only a narrow smoke check; at
full scope they are 2.5% of the total. Simulation is most of it, and
*initialisation is more than a quarter* -- almost all of that the `large`
system, which spends 54 s on setup before its first MD step, and pays it again
on every repeat.

That is because initialisation scales far worse than simulation does:

| System | Atoms | Initialisation | Per MD step |
|---|---|---|---|
| `tiny` | 11,904 | 0.34 s | 8.2 ms |
| `medium` | 45,530 | 4.24 s | 33.9 ms |
| `large` | 172,541 | 54.05 s | 136.8 ms |

14.5x the atoms costs 15.4x per step -- essentially linear -- but **159x the
initialisation**. `track_init_walltime` exists to keep an eye on it.

So the levers, in order of value:

1. **Drop the `large` tier** for routine runs. Half the cost of a one-variant
   sweep, and most of that is setup you are not measuring. Run it weekly.
2. **Use `quick.conf.json`** -- one variant instead of four.
3. **Fewer repeats.** 3 -> 2 costs some robustness against an unlucky run; the
   metric is the minimum, so it degrades gracefully.
4. **`--skip-existing-successful`** so repeat invocations do not redo work.

What *not* to do:

* **Do not shorten NSTLIM to save time.** Measured spread over five runs of the
  `tiny` system: 0.55% at 4000 steps, 1.13% at 250, 3.73% at 125. Longer runs
  average out noise, which is the whole point of the step count. Shorten the
  scope, never the measurement.
* **Do not use `asv --parallel`.** Concurrent runs contend for cores. A single
  build running alongside a benchmark inflated one measurement here by 11%, and
  another by 54%.

```sh
# fast smoke check
GROMOS_BENCH_TIERS=tiny GROMOS_BENCH_THREADS=4 GROMOS_BENCH_REPEATS=1 \
  ./asv run --config quick.conf.json HEAD^!

# routine sweep: coarse first, then fill in -- --skip-existing-successful
# means the staged version costs no more than running the last pass alone,
# but leaves a usable timeline at every stage
export GROMOS_BENCH_TIERS="tiny medium" GROMOS_BENCH_REPEATS=2
for s in 17 34 67 134; do
  ./asv run --config quick.conf.json --skip-existing-successful \
            --steps $s 023d27db..master
done
```

### Paths

Everything asv writes lives under `.asv/` beside the config (gitignored), so the
suite runs on any machine without editing:

```
.asv/env       virtualenvs, project checkouts, cached builds
.asv/results   the results database, one directory per machine
.asv/html      the generated site
```

Expect `.asv/env` to reach several GB: a cached build is a ~55 MB binary,
`build_cache_size` keeps 8 of them, and there is one environment per variant.

Two paths are *not* in the config, because they must point at local disk that
suits the machine rather than the repository:

| Variable | Default | Purpose |
|---|---|---|
| `GROMOS_BENCH_SYSTEMS` | `/local/gromos/bench_systems` | staged inputs (~39 MB) |
| `GROMOS_BENCH_SCRATCH` | `/local/gromos/bench_scratch` | per-run working files |

Set both on a machine without `/local/gromos`. They are deliberately kept out of
the asv matrix so they do not become part of every environment's name and hash.

### Running on another machine

1. Clone the repository. Everything code-related is committed.
2. Install the system dependencies (see [Prerequisites](#prerequisites)).
3. Install asv:
   `python3 -m pip install --target /local/gromos/bench_tools asv virtualenv`
   (set `GROMOS_ASV_TOOLS` if you put it elsewhere).
4. Provide the benchmark systems, either by staging them from the upstream set
   with `python3 tools/stage_systems.py`, or by copying the cache across and
   running `python3 tools/stage_systems.py --verify`. Verify either way:
   identical inputs are what make two machines' results comparable.
5. `./asv machine --yes`, then run as usual.

Results from several machines can share one `results/` directory -- they are
keyed by machine name and published as selectable series. `env/` is
machine-local and must never be shared.

Note that the compiler is part of the measurement. A second machine on a
different gcc produces a series that is internally consistent but not
comparable with this one; use `GROMOS_CXX` to force a match if you intend to
read across them.

### Prerequisites

* `nvcc` for the `cuda` variant. It is frequently not on `PATH`; set
  `GROMOS_NVCC` to its location.
* FFTW, from the system packages. Note that the build deliberately uses the
  system compiler even when a conda environment is active: a conda toolchain
  does not see the system FFTW, and more importantly the compiler is part of
  what the timeline measures. Set `GROMOS_CXX` to override.
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
| `GROMOS_CXX` | C++ compiler (default: the system one, not conda's) |
| `GROMOS_CUDA_LIB` | dir holding `libcudart` (default: `<nvcc>/../../lib`) |

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
quick.conf.json       the same, restricted to one variant, for timelines
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
