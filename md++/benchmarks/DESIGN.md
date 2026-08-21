# How the benchmark suite works

`README.md` covers how to *use* the suite. This describes how it is built and
why, for anyone extending it or debugging a run that goes wrong.

The short version: asv is a Python benchmarking tool, md++ is a C++ program, and
most of the work here is bridging that gap. asv contributes the commit-walking,
the results database and the web UI; everything about actually building and
measuring md++ is supplied by this directory.

## Contents

- [The central idea](#the-central-idea)
- [Layers](#layers)
- [What happens during a run](#what-happens-during-a-run)
- [Build integration](#build-integration)
- [Configuring a simulation](#configuring-a-simulation)
- [Reading the results](#reading-the-results)
- [The parameter grid](#the-parameter-grid)
- [Measurement discipline](#measurement-discipline)
- [Guard rails](#guard-rails)
- [Extending the suite](#extending-the-suite)
- [asv behaviours worth knowing](#asv-behaviours-worth-knowing)

## The central idea

**md++ measures itself; the harness only reads what it reports.**

At the end of every run, `program/md.cc` prints a timing summary, and
`md.print_timing()` emits a `TIMING` block produced by
`util::Algorithm_Timer::print` (`src/util/timing.cc:377`) that breaks the step
down per algorithm, with named sub-timers underneath the ones that have them:

```
TIMING
          NonBonded                              0.031     100.00%
             - pairlist                                    0.009  28.22%
             - shortrange solvent-solvent                  0.008  27.58%
          Shake                                  0.009     100.00%
END

Wall time initialisation (s):       0.660
Wall time simulation (s):           1.619
Wall time total (s):                2.279
Performance (ns/day):              21.341
```

Two consequences follow, and they shape everything else.

**Each configuration runs exactly once, and many metrics are scraped from that
one run.** So the benchmarks are asv `track_*` functions, not `time_*` ones.
asv's timing machinery is designed for microbenchmarks: it would re-run each
configuration repeatedly and measure the whole process, so a multi-minute MD run
would be timed together with process startup, topology parsing and the input
coordinate read. Running once and reading md++'s own numbers is both cheaper and
more precise.

**The per-algorithm series are the point.** A single wall-clock number tells you
a commit got slower. The breakdown tells you the pairlist got slower while the
short-range loop did not, which is the difference between a graph you look at
and one you can act on.

## Layers

```
   asv  ──────────────────  commit walking, results database, web UI
    │
    │  build_command / install_command          (asv.conf.json)
    ▼
   scripts/asv_build.sh ───  cmake configure / build / install
   scripts/asv_install.sh    install prefix -> asv environment
    │
    │  benchmark_dir
    ▼
   benchmarks/md_bench.py ─  PlainMD: the asv-facing benchmark class
    │
    ├── _systems.py ───────  which system, which files, how many steps
    ├── _imd.py ───────────  rewrite the simulation input
    ├── _runner.py ────────  launch md++, per variant, pinned
    └── _omd.py ───────────  parse md++'s output into metrics
```

Each module below `md_bench.py` is independently usable, which is what
`tools/calibrate.py` relies on to exercise the harness without asv involved.

| File | Lines | Role |
|---|---|---|
| `asv.conf.json` | 64 | asv configuration; the four-variant build matrix |
| `nightly.conf.json` | 59 | the same, restricted to one variant, for timelines |
| `asv` | 44 | wrapper supplying `PYTHONPATH` and the working directory |
| `scripts/asv_build.sh` | 161 | cmake build, variant flags, build-time guards |
| `scripts/asv_install.sh` | 26 | copy the install prefix into an asv environment |
| `benchmarks/md_bench.py` | 188 | the `PlainMD` benchmark class |
| `benchmarks/_runner.py` | 226 | variant dispatch, pinning, repeats |
| `benchmarks/_omd.py` | 126 | omd parser |
| `benchmarks/_imd.py` | 135 | block-aware imd reader/patcher |
| `benchmarks/_systems.py` | 120 | system registry and NSTLIM table |
| `tools/stage_systems.py` | 140 | stage inputs into the local cache |
| `tools/calibrate.py` | 92 | measure per-step cost; standalone smoke test |

## What happens during a run

For each commit and each build variant:

1. **asv checks the commit out** into `<env>/project`.
2. **asv calls `scripts/asv_build.sh`**, which runs cmake and installs into
   asv's per-commit build cache.
3. **asv calls `scripts/asv_install.sh`**, which copies that prefix into the
   environment root, putting the binary at `<env>/bin/md`.
4. **asv imports `md_bench.py`** and calls `setup_cache()` once.
5. `setup_cache()` runs every selected (system, threads) combination:
   `_systems` resolves the input files, `_imd` writes a benchmark-specific imd,
   `_runner` launches md++, `_omd` parses the output. Results accumulate in a
   dict keyed `"<tier>/<threads>"`.
6. **asv calls each `track_*` function**, passing that dict as the first
   argument. Each is a dictionary lookup — no simulation runs here.
7. **asv writes** one JSON file per commit per environment under `results/`.

Step 4 is the crux. asv calls `setup_cache` once per benchmark class and hands
its return value to every benchmark method, so all eight metrics come from the
same set of runs. Doing the work in `setup()` instead would re-run every
simulation once per metric — eight times the cost, and eight mutually
inconsistent sets of numbers.

The returned dict must be picklable; it holds only floats.

## Build integration

asv's defaults (`python -m build`, `pip install <wheel>`) are meaningless for a
C++ program, so both are replaced:

```json
"build_command":   ["{conf_dir}/scripts/asv_build.sh {build_dir} {build_cache_dir}"],
"install_command": ["{conf_dir}/scripts/asv_install.sh {build_cache_dir} {env_dir}"],
"uninstall_command": []
```

Three decisions are embedded here.

**Why a script rather than cmake directly in the JSON.** asv runs build commands
*without a shell*: it `shlex.split`s the string and substitutes only
`{placeholder}` tokens (`asv/util.py`, `asv/environment.py`). `$GROMOS_VARIANT`
would never expand. The variant reaches the script through its environment
instead, which asv does pass through (`extra_env=self.build_env_vars`).

**Why `{conf_dir}` and not `{build_dir}`.** `{build_dir}` is the checked-out
historical commit, which does not contain this benchmark suite — every commit
older than the suite would fail with "no such file". `{conf_dir}` is your
working checkout, so the hooks always exist. This is what makes benchmarking
pre-existing history possible at all.

**Why the variant flags live in the script.** `scripts/asv_build.sh` maps
`GROMOS_VARIANT` to cmake flags:

| Variant | cmake flags | Binary | Launch |
|---|---|---|---|
| `serial` | *(none)* | `md` | direct, 1 thread |
| `omp` | `-DOMP=ON` | `md` | direct, `OMP_NUM_THREADS` |
| `cuda` | `-DOMP=ON -DCUKERNEL=ON` | `md` | direct + `INNERLOOP` CUDA |
| `mpi` | `-DMPI=ON` | `md_mpi` | `mpirun -np N` |

Deriving them there keeps the asv matrix one-dimensional and makes the
combinations md++ rejects (`cmake/options.cmake`: OMP with MPI, CUKERNEL without
OMP) impossible to express.

Note that `cuda` is an OpenMP build. Whether the CUDA kernel is *used* is a
runtime choice made in the imd, not a build-time one — see below.

### Installing a C++ build into a Python environment

`scripts/asv_install.sh` is `cp -a <prefix>/. <env_dir>/`. Because an asv
virtualenv's root is `sys.prefix`, the binary lands at `sys.prefix/bin/md`, and
`_runner.find_binary()` finds it with no path configuration and no knowledge of
which environment it is in. Each variant automatically resolves to its own
build.

### Pinning the toolchain

`scripts/asv_build.sh` sets `CMAKE_CXX_COMPILER` explicitly rather than
inheriting whatever the invoking shell offers. The visible reason is that an
activated conda environment puts its own toolchain first, and that toolchain
does not see the system FFTW, so configuration fails with `fftw3 library could
not be found`.

The real reason is the case where it *does not* fail. The compiler is part of
the measurement: a timeline built partly with the system gcc and partly with a
conda gcc compares compilers rather than commits, and nothing in the recorded
results would reveal it. `GROMOS_CXX` overrides the choice deliberately, which
is the supported way to compare compilers on purpose.

Two things must survive that copy:

- **RPATH.** cmake gives the build tree an RPATH covering what it linked and
  *strips it on install*. With the CUDA runtime living in a conda environment
  rather than on the loader path, the installed binary then dies with
  `error while loading shared libraries: libcudart.so.12` even though the same
  binary runs fine in the build tree. `CMAKE_INSTALL_RPATH_USE_LINK_PATH` is not
  enough — it only covers libraries the link line names by full path, and the
  CUDA runtime arrives via a `-L` search path — so `CMAKE_INSTALL_RPATH` is set
  explicitly from `GROMOS_NVCC`.
- **Optimisation.** See [Guard rails](#guard-rails).

## Configuring a simulation

Benchmark inputs are three equilibrated protein-in-solvent systems staged into a
local cache by `tools/stage_systems.py`. They are not committed: they total
~39 MB, and a manifest of sha256s is verified before every run, because inputs
that drift silently turn the timeline into a measurement of nothing.

`_imd.prepare_benchmark_imd()` derives a benchmark input from the system's
production imd, changing exactly three blocks:

| Block | Change | Why |
|---|---|---|
| `STEP` | NSTLIM from the registry, start time 0 | target a run duration; make runs comparable |
| `WRITETRAJ` | all frequencies to 0 | otherwise the benchmark measures the filesystem |
| `INNERLOOP` | `0 0` for CPU, `4 0 1 0` for CUDA | selects the nonbonded kernel at *runtime* |

The `INNERLOOP` line is why the parser is block-aware rather than line-based.
Its parameter count *varies with the method*: with `NTILM=0` only two values may
appear, and supplying the four a CUDA run needs makes md++ fail during
initialisation with

```
ERROR In_Parameter : Left-over parameters in INNERLOOP block
```

A naive substitution preserving the original column count silently produces an
unrunnable input for CPU runs. `_imd.ImdFile` therefore parses the file into
ordered `NAME -> [lines]` blocks, preserves the `#` comment lines that document
each block's columns, and lets `set_values()` change the parameter count freely.

Because the kernel is selected here, one CUDA-capable build serves both CPU and
GPU benchmarks.

## Reading the results

`_omd.parse()` extracts, from md++'s own output:

| Key | Source |
|---|---|
| `wall_initialisation`, `wall_simulation`, `wall_total` | summary lines |
| `ns_per_day` | `Performance (ns/day)` |
| `timer.<Algorithm>` | `TIMING` block, main timers |
| `subtimer.<Algorithm>.<name>` | `TIMING` block, indented sub-timers |

Sub-timer names contain spaces (`longrange solvent-solvent`), so the name is
matched non-greedily up to the numeric columns rather than by splitting on
whitespace. `unaccounted time` is skipped — it is bookkeeping, not an algorithm.

The exported metrics, and two deliberate choices:

- **`track_sim_walltime` uses `Wall time simulation`, never `total`.** For the
  large system, reading the 29 MB input coordinates takes about four times
  longer than the simulation being measured, so total wall time would be
  dominated by I/O unrelated to the integrator.
- **`track_ns_per_day` is informational only.** asv's regression detection
  assumes lower-is-better, so on a higher-is-better series every genuine speedup
  would be reported as a regression.

## The parameter grid

`PlainMD.params` is **fixed**:

```python
params = [list(TIERS), list(THREADS)]      # [tiny|medium|large] x [1,2,4,8,16]
```

It must never depend on the environment. asv keys results by *position* in this
grid and records the grid in `results/benchmarks.json`; if it changed shape
between invocations, results written under one grid would be orphaned by the
next run under another, and `asv publish` would emit graphs containing literally
nothing while reporting success.

`GROMOS_BENCH_TIERS` and `GROMOS_BENCH_THREADS` therefore narrow only what
*runs*. Combinations left out are absent from the `setup_cache` dict, and
`setup()` reports them as skipped. A quick `GROMOS_BENCH_TIERS=tiny` check and a
full nightly run write mutually compatible results.

Skipping is done by raising `NotImplementedError` **from `setup()`**, not from
the benchmark function. asv catches it only in `do_setup`
(`asv_runner/benchmarks/_base.py`); the same exception raised inside a benchmark
method is recorded as a *failure*. That is also how the serial variant declines
every thread count above one.

## Measurement discipline

MD benchmarking is easy to get wrong, so several things are deliberate:

- **No trajectory writing** (`WRITETRAJ` zeroed) and a scratch directory on
  local disk, so the filesystem is not in the measurement.
- **Initialisation excluded** — `Wall time simulation`, never `total`.
- **Threads bound to cores** (`OMP_PROC_BIND=close`, `OMP_PLACES=cores`) and the
  process pinned to a fixed CPU set with `taskset`, so a run is not migrated
  mid-measurement. Disable with `GROMOS_BENCH_PIN=0`.
- **Minimum across repeats, not mean.** Interference from the rest of the
  machine can only ever make a run slower, so the fastest observation is the one
  closest to the machine's true capability. All metrics come from the single
  fastest repeat rather than being reduced independently, which keeps the
  per-algorithm breakdown internally consistent.
- **NSTLIM sized per system** so a run lands near 30 s at four threads.
  `tools/calibrate.py` re-measures and suggests values after a hardware change.

Two things the suite cannot fix and only warns about: the CPU frequency governor
(set it to `performance`; `calibrate.py` warns otherwise), and other load on the
machine. With a single repeat on a `powersave` governor the run-to-run spread is
around 1.4%, which is the floor on what a regression must exceed to be visible.

`GROMOS_BENCH_TRACE=<file>` appends one line per simulation recording pid,
variant, thread count, effective `OMP_NUM_THREADS`, CPU affinity and timing.
It is the quickest way to confirm what actually ran.

## Guard rails

Failures that would otherwise be recorded as plausible-looking numbers:

| Guard | Where | Catches |
|---|---|---|
| Debug build rejected | `scripts/asv_build.sh` re-reads `CMAKE_BUILD_TYPE` from `CMakeCache.txt` after building | `GROMOS_CMAKE_ARGS` overriding the build type |
| Debug build rejected | `_omd.parse()` looks for md++'s own "This is a debug code" banner | any debug binary, however produced |
| Pre-rename commit | `scripts/asv_build.sh` checks for `md++/` vs `gromosXX/` | commits older than `284700de` (2023-10-21) |
| MPI variant, non-MPI binary | `_runner.find_binary()` | `mpirun` launching N *independent* simulations that compete for cores and look like a severe slowdown rather than an error |
| Missing / altered inputs | `_systems.load_manifest()` at the start of `setup_cache` | a drifting timeline |
| Unpinned compiler | `scripts/asv_build.sh` sets `CMAKE_CXX_COMPILER` | an activated conda toolchain silently replacing the system one |
| Run failed | `_omd.parse()` requires `MD++ finished successfully` | any failure, reported with the command, exit code and stderr tail |

A debug build is several times slower, so without those two guards it would
appear on the timeline as an enormous regression rather than as a mistake. The
`_omd` check is the stronger of the pair: it holds regardless of how the binary
was produced, including one supplied via `GROMOS_MD_BINARY`.

## Extending the suite

**Another metric from an existing run.** Add a `track_*` method to `PlainMD`
reading another key from the cache dict — no extra simulations. Give it a
`.unit`. Note that changing what an existing metric measures makes old and new
points on the same graph incomparable; `./asv rm` that series.

**Another system.** Add an entry to `_systems.SYSTEMS` with `pdb`, `atoms` and
`nstlim`, then stage it. Adding a tier changes the parameter grid, which
orphans existing results for the reason given above.

**Another method (AEDS, replica exchange, ...).** Add a benchmark class beside
`PlainMD`. If it needs different imd blocks, extend `_imd` — `ImdFile` is
general. If it needs a different binary, extend `_runner.binary_name()`.

**Another build variant.** Add it to `VARIANTS` in `_runner.py`, to the `case`
in `scripts/asv_build.sh`, and to the matrix in `asv.conf.json`.

**Per-variant step counts.** A system may carry
`"nstlim_by_variant": {"serial": 1500}`, which `_runner` consults before falling
back to `nstlim`. Changing either value breaks comparability with results
already recorded.

## asv behaviours worth knowing

Things that cost real debugging time here, all confirmed against the installed
asv source:

- **`asv publish` only graphs commits reachable from a branch listed in
  `branches`.** Everything else is dropped *silently* — it logs
  `Couldn't find <hash> in branches (...)` and then reports success having
  generated nothing. This is the usual cause of "No graphs to load".
  `asv show` and `asv compare` ignore `branches` and read the results directory
  directly, so they still work for unlisted branches.
- **asv chdirs to the directory containing the config file** (`asv/main.py`) and
  resolves `repo`, `benchmark_dir`, `env_dir`, `results_dir`, `html_dir` and
  `{conf_dir}` relative to *that*, not to where you invoked it. This is what lets
  every path in the config be relative, and why `nightly.conf.json` must sit
  beside `asv.conf.json`. A config kept elsewhere needs absolute paths
  throughout.
- **`-E/--environment` selects only the environment type and Python version.**
  It cannot filter matrix variables, so restricting to one build variant
  requires a separate config file.
- **Build commands get no shell** — `shlex.split` plus `{placeholder}`
  substitution only. Environment variables reach them through the process
  environment, not through `$VAR` in the JSON.
- **`NotImplementedError` is a skip signal only from `setup`.** From a benchmark
  method it is a failure.
- **The benchmark code always comes from your working checkout**, never from the
  commit under test. Improving the benchmarks and re-running rewrites history
  with the new code; benchmark changes are not versioned alongside results.
- **A failed build leaves an empty `prefix/bin`** in the build cache, which is a
  quick way to tell a build failure from a benchmark failure.
- asv eagerly loads its mamba plugin and logs an error when `libmambapy` is
  absent; its own guard misses the module because it is named `_mamba_helpers`
  rather than `mamba`. Harmless, and filtered by the `./asv` wrapper.
