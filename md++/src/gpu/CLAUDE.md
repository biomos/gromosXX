# CLAUDE.md — GROMOS CUDA work

This file guides work on GROMOS's CUDA/GPU acceleration layer. It's scoped to
`src/gpu` but touches several places outside it (see "Cross-cutting files"
below), because GPU support is being woven into the core simulation classes'
call paths, not bolted on beside them.

Read `src/gpu/PLAN.md` before making architectural changes — it's the
decision record for this effort and takes precedence over anything in
`src/gpu/notes.md` (which is left as historical brainstorming; still useful
for mechanical how-tos, but not for decisions it left open or that PLAN.md
has since settled differently).

`src/gpu_legacy` is the previous CUDA implementation, being replaced. Don't
extend it. It does contain real, working prior art for the pairlist/cell-list
problem (`interaction/nonbonded/pairlist/grid_pairlist.cu`) with measured
performance characteristics in its comments — worth reading before
reimplementing that piece, not worth reusing verbatim (different container
types).

## The one-sentence philosophy

Minimize host↔device transfers, not per-kernel microseconds. The performance
goal is running as much of one MD step on GPU as possible so data never
leaves the device mid-step — that's why several decisions below (GPU mirror
caching, the hard error on CPU-pairlist-under-CUDA, the fallback warning)
exist: they're all about making transfer cost visible and avoidable, not
about squeezing individual kernels.

## Data residency rule: GPU stays authoritative

This is a hard rule, not a per-algorithm judgment call — established after a
session of finding and fixing violations one at a time (see git log around
`cuda_shake`/`cuda_settle`/`constraint_force_publish_kernels`/
`virial_accumulate_kernels` for the concrete before/after).

**Any simulation data (positions, velocities, force, constraint_force,
virial_tensor, energies, box, ...) that a GPU-native algorithm produces or
consumes must live in the shared `CudaManager` mirror
(`sim.cuda().configuration_view()`/`mark_gpu_dirty()`), and stay there.**

- If multiple GPU-native algorithms accumulate into the same shared
  quantity (virial_tensor, force), the accumulation itself happens on-device
  — `atomicAdd` into the mirror's own field, not each algorithm's private
  buffer downloaded and summed on the host. `gpu/cuda/memory/
  virial_accumulate_kernels.h` is the reference pattern: one small
  `atomicAdd`-9-doubles kernel, called by every contributor, into a field
  zeroed once per step (`CudaManager::zero_mirror_force()`).
- If the writers' target ranges are disjoint (constraint_force: solute vs
  solvent constraint algorithms never overlap), a plain on-device write is
  enough — no atomics needed, but it's still a device-to-device operation,
  never a host round trip. `gpu/cuda/algorithm/constraints/
  constraint_force_publish_kernels.h` is that pattern.
- A GPU-native algorithm may still keep a **private** per-term GPU buffer
  for its own internal use (e.g. each bonded term's own energy/virial
  before merging into the shared total) — that's fine, even expected, when
  the decomposition itself is meaningful (per-interaction-type energies are
  reported separately in output). What's not fine is downloading that
  private buffer to CPU and merging it there every step "just because it's
  easier" — the merge must happen on-device too, into the shared mirror
  field, via the on-device patterns above.
- **The CPU never gets an eager, unconditional publish.** A CPU-backend
  algorithm that needs to read GPU-resident data is responsible for that
  publish itself, through the exact same mechanism every other consumer
  uses: either its own `gpu_mirror_touches()` mask (if it's registered in
  `Algorithm_Sequence`, e.g. `Pressure_Calculation`'s
  `gpu_mirror_touches() == MIRROR_VIRIAL`, whose before-hook `flush_gpu_
  dirty()` is the *only* thing that copies virial_tensor to the host, once
  per step, only because that one CPU-only reader needs it) or an explicit
  `sim.cuda().flush_gpu_dirty(conf, <fields>)` call at the point of actual
  need (trajectory/checkpoint output, `program/md.cc`'s output-cadence
  flush is the reference example). Never flush "just in case" on a fixed
  schedule unrelated to an actual reader.
- This includes test code: a test that calls `apply()`/
  `calculate_interactions()` standalone (bypassing `Algorithm_Sequence::
  run()`'s automatic per-algorithm flush) and then reads `conf` directly
  must call `flush_gpu_dirty()` itself before that read, mirroring what a
  real consumer would do — this is not a workaround, it's the test
  correctly exercising the same contract a real caller has to honor. Don't
  "fix" a stale-read test failure by making the algorithm publish eagerly
  instead.
- Before writing a CPU-merge loop (`conf.old().foo(i) += ...`) or a plain
  host readback for anything GPU-computed, stop and ask whether this
  belongs in the mirror instead — the default assumption for new
  GPU-native code should be "stays on GPU," not "syncs back are cheap
  enough."

## Current status (see PLAN.md §2 for the full decision log)

- Backend dispatch (`util::cpuBackend`/`gpuBackend`, `algorithm::AlgorithmB<Backend>`,
  `algorithm::make_algorithm<AlgT>`) is implemented and works — validated by
  `Remove_COM_Motion`, `Lattice_Shift_Tracker`, `Leap_Frog_Velocity/Position`.
  Use this pattern for every new GPU-portable algorithm.
- `gpu::cuvector`/`cuhvector`, `gpu::Container`, `gpu::Interaction_TileT`/
  `TileVecT`/`TileContainerT`, `gpu::Periodicity<BOUNDARY>` are solid
  infrastructure — build on them, don't replace them.
- `CudaManager`'s worker-thread/task-queue layer and its "shallow copy with a
  printed warning" semantics are being removed (PLAN.md §3.1) — don't copy
  that pattern into new code; new CUDA-resource-owning classes are move-only.
- `gpu::Topology`/`gpu::Configuration` are being moved out of
  `topology::Topology`/`configuration::Configuration` into an identity-keyed
  cache inside `CudaManager`, reached via `sim.cuda()` (PLAN.md §3.2). If you
  see `conf.copy_to_gpu()`, `conf.m_gpu`, `topo.get_gpu_view()` in code you're
  touching, it's using the old (being-removed) API — migrate it, don't
  extend it.
- `CudaMemoryManager` has several known half-finished spots (uninitialized
  `m_device_id`, `allocate<T>(count)` ignoring `count`, a smart-pointer type
  that's never actually wired to the manager) — see PLAN.md §3.3 for the
  checklist. Fix these before building more on top.
- The pairlist is the big remaining piece of real work: exactly one GPU-
  native tile/cell-list implementation, no user-facing pairlist choice under
  CUDA (PLAN.md §6). The commit(s) after `a5118b63d` that build a CPU
  pairlist + per-step full upload + one-thread-per-pair kernel are being
  discarded (PLAN.md §7) — don't extend `cuda_nonbonded_interaction.cc` or
  `nb_kernels.cu` as they stand; their RF-math bodies are salvageable but the
  data flow around them is not.

## Cross-cutting files (outside `src/gpu`)

- `simulation/simulation.h` — owns `m_cuda_manager` unconditionally; `sim.cuda()`
  is the single access point for all things CUDA, including the GPU mirrors
  once §3.2 lands.
- `topology/topology.h`, `configuration/configuration.h` — currently include
  `gpu/cuda/...` headers and hold GPU-mirror members; this coupling is being
  removed (PLAN.md §3.2). Don't add more GPU awareness to these classes.
- `algorithm/algorithm.h`, `util/backend.h` — the backend-dispatch machinery.
  `make_algorithm<AlgT>`'s CPU-fallback-when-no-GPU-backend path needs a
  warning added (PLAN.md §8) — check whether that's landed before assuming
  silent fallback is intentional.
- `interaction/nonbonded/pairlist/*` — six CPU `Pairlist_Algorithm`
  subclasses exist (`Standard_*`, `Standard_*_Atomic`, `Extended_Grid_*`,
  `Grid_Cell_*`), all producing the same `PairlistContainer` shape. The GPU
  one is a seventh, selected exclusively when `accelerator == gpu_cuda`
  (never user-choosable alongside the CPU ones).
- `simulation/parameter.h` — `plist_struct` (the `PAIRLIST` block) gets a new
  `skin` field (PLAN.md §6.3); `gpu_struct` (accelerator selection) is
  unchanged.

## Conventions to follow

- **CPU/CUDA split**: `.cc` files compile in CPU-only builds and use the
  `DISABLED()`/`DISABLED_VOID()`/`DISABLED_MSG()` macros (`cuheader.h`) for
  every method body. `.cu` files compile only under `USE_CUDA` and hold the
  real implementation. Never put `#ifdef USE_CUDA` blocks inside a single
  shared source file if a `.cc`/`.cu` split is possible — that's the whole
  point of the split (see `notes.md`'s "Code splitting for CPU and GPU").
- **New GPU-portable algorithm**: turn the class into
  `template <typename Backend = util::cpuBackend> class Foo : public Algorithm, private AlgorithmB<Backend>`,
  add `is_supported_backend`, split the implementation into `foo_cpu.cc` +
  `foo_gpu.cu` (or use `if constexpr` in one `.tcc` if the difference is
  small). Instantiate via `algorithm::make_algorithm<Foo>(sim, ...)`, never
  by naming `Foo<util::gpuBackend>` directly in code that must compile
  without `USE_CUDA`.
- **Precision**: use `FPL_TYPE`/`FPH_TYPE` (and `FPL2/3/4/9_TYPE`,
  `FPH2/3/4/9_TYPE`) from `gpu/cuda/memory/precision.h`, never bare
  `float`/`double`, in any code that's meant to work across all three
  `FP_PRECISION` settings. Default build is mixed precision (`FP_PRECISION=2`)
  — code should be correct at `=1` and `=3` too since those remain supported
  debug/comparison builds.
- **Memory ownership**: algorithm-private persistent GPU data is an ordinary
  class member (`gpu::cuvector<T>`, or `CudaMemoryManager::cumm_unique_ptr<T>`)
  — RAII, no central registry call needed. Only go through `CudaManager`/
  `CudaMemoryManager` for things that are genuinely shared (device selection,
  the Topology/Configuration mirrors, the pinned staging pool).
- **Boundary conditions**: dispatch on `math::boundary_enum` via
  `gpu::Periodicity<BOUNDARY>` template parameter + `SPLIT_BOUNDARY`-style
  runtime switch at the host call site, matching how the CPU code already
  does it (see `math::boundary_implementation.h` usage patterns) — don't
  invent a second boundary-dispatch mechanism.
- **No copyable CUDA-resource owners.** Any class holding a `cudaStream_t`,
  a device pointer it's responsible for freeing, or a worker/thread must be
  move-only (`= delete` the copy ctor/assignment). Don't add a "shallow copy
  + printed warning" pattern anywhere, even under time pressure — it's a
  double-free waiting to happen, not a style shortcut.

## Testing expectations

Every GPU code path needs a correctness test against the CPU reference before
it's considered done — see `PLAN.md` §9. The pairlist specifically needs an
equivalence test (same pair set, including exclusions, as
`Standard_Pairlist_Algorithm`) before its output is trusted for any force/
energy comparison. Don't hand-wave this with a manual "looks about right"
check — plug it into the existing GROMOS regression-test infrastructure.

## What not to do

- Don't reintroduce a `#ifdef USE_CUDA` scattered through shared logic when a
  `.cc`/`.cu` split or a `Backend` template parameter would do instead.
- Don't add a new pairlist-algorithm choice under CUDA. There is exactly one.
- Don't build a "CPU pairlist then copy to GPU" fallback as a permanent,
  always-available code path — that's the discarded approach, and if it's
  always available it becomes the silent default under any misconfiguration.
- Don't add worker threads / task queues to the CUDA manager layer. Use
  streams and events.
- Don't put GPU-mirror state back into `topology::Topology` or
  `configuration::Configuration`.

## Git workflow

- Always work on a branch named `cuda_claude_<short_feature_name>`, never commit directly to `master`
- Decide on checkpoints and make separate branch for every separately reviewable part
- The branch at its final stage should compile properly and pass all the tests
- Commit after each logical, working change with a clear message
- Never force-push
- Never merge to `master` yourself — leave that for human review
- Tests are not required to pass on every commit, but the build must always succeed. The full test suite must pass before the branch is considered complete (see Standard workflow).

## Build

This project uses **CMake with the Ninja generator**, driven through CMake
presets — not autotools/make, despite what older notes or muscle memory
might suggest. There is no `configure` script and no `Makefile` in the
build directory; `make` will fail with "No targets specified and no
makefile found." Use `ninja` or `cmake --build`, never `make`, in this repo.

​```bash
cmake --preset cuda-on
cmake --build --preset cuda-on -j $(nproc)
ctest --preset cuda-on --output-on-failure
​```

If a build preset isn't defined for a given configure preset, build directly
from the build directory instead:
​```bash
cd build
ninja -j $(nproc)
ctest --output-on-failure
​```

Rules:
- Always check `cmake --list-presets` first if unsure which preset to use;
  presets are self-documenting via their `displayName`/`description` fields
  — read those rather than guessing which one fits the task.
- Do not re-run `cmake --preset ...` on every build — only if `CMakeLists.txt`,
  a `cmake/*.cmake` module, or the preset file itself changed. Plain `ninja`
  (or `cmake --build`) is sufficient for a normal edit/rebuild cycle.
- If you changed configure-time options (e.g. toggling `USE_CUDA`), delete
  the build directory and reconfigure from scratch rather than reusing a
  stale cache: `rm -rf build && cmake --preset cuda-on`. CMake does not
  reliably overwrite a cached variable just because the preset changed.
- A build is only "done" when `cmake --build ...` exits 0 AND
  `ctest ... --output-on-failure` exits 0. Do not report success otherwise.

## Build/test failure procedure

If `cmake --preset ...` (configure step) fails:
1. Read the exact error in the output (not just the final summary line).
2. Check `build/CMakeFiles/CMakeConfigureLog.yaml` (CMake ≥3.26) for the
   specific failing check — this is the CMake equivalent of autotools'
   `config.log`. On older CMake, `build/CMakeCache.txt` shows what was
   actually detected/set, though not the full narrative log.
3. Do not disable a feature (e.g. drop `USE_CUDA`) to work around a missing
   dependency — report what's missing instead, unless explicitly told the
   dependency is optional.
4. `CMAKE_CUDA_ARCHITECTURES must be non-empty if set` means the variable is
   set-but-empty, not merely unset — check `CMakeLists.txt` for where it's
   set relative to `enable_language(CUDA)`. It must be set *before* that
   call, not after; CMake validates/locks it in during `enable_language`.

If the build step (`ninja` / `cmake --build`) fails:
1. Fix the actual compile/link error shown — do not suppress warnings-as-errors
   by weakening flags unless the warning is a false positive you can justify
   in the commit message.
2. Re-run the build after each fix; do not batch multiple unverified fixes
   before re-testing.

If `ctest` fails:
1. Identify which specific test(s) failed: `ctest --output-on-failure` shows
   the failing test's output directly; `build/Testing/Temporary/LastTest.log`
   has the full detail if you need more.
2. Distinguish a real bug (code produces wrong result) from a broken test
   (test itself is outdated/incorrect) — fix the actual bug by default; only
   modify a test if you can justify the test itself is wrong.
3. Re-run only the failed test first if possible
   (`ctest -R <test_name> --output-on-failure`), then the full suite, before
   considering it fixed.
4. If a fix is non-obvious after 2-3 attempts, stop and leave a clear note in
   the commit/PR description rather than making increasingly speculative
   changes.

## Writing new tests

- New tests go in `src/check`
- Follow the existing test file naming/structure — look at 2-3 existing tests
  before writing a new one
- Register new test files in the relevant `CMakeLists.txt`
  (`add_test(...)` / the test list for that directory) so `ctest` picks them
  up automatically
- Every new feature should come with an appropriate test
- Every bug fix should come with a regression test that fails before the fix
  and passes after, where practical
- Run `ctest --preset cuda-on --output-on-failure` (or `ctest` from `build/`)
  after adding a test to confirm it's actually being picked up by the
  harness (not silently skipped) — check `ctest -N` to list discovered
  tests without running them, if you just want to confirm registration

## CUDA build notes
- Verify `nvcc --version` succeeds before building; if CUDA isn't found, do
  not silently drop `USE_CUDA`/`--with-cuda` — report it.
- `CMAKE_CUDA_ARCHITECTURES` is set explicitly in `CMakeLists.txt`
  (currently `80 86 89 90 120` — `120` added for consumer Blackwell cards,
  e.g. RTX 50-series, compute capability 12.0, after a real "PTX was
  compiled with an unsupported toolchain" failure running the binary on a
  driver that couldn't JIT-compile PTX for that architecture; baking in
  native SASS avoids depending on the runtime driver's JIT compiler at
  all) — do not rely on preset-level overrides for this; the in-source
  `set()` runs before `enable_language(CUDA)` and takes effect regardless
  of what a preset provides. If you need to target a different/additional
  architecture, edit it there, not in a preset.

## Standard workflow
1. Make code change
2. Build (`cmake --build --preset cuda-on -j $(nproc)`, or `ninja` from `build/`)
3. Test (`ctest --preset cuda-on --output-on-failure`, or `ctest` from `build/`)
4. If failures: diagnose → fix → rebuild → retest (do not skip straight to committing)
5. A commit requires a successful build, but not a fully passing test suite - WIP commits on a feature branch may carry known, understood test failures.
6. Every test failure must still go through the Build/test failure procedure
   (diagnose before deferring) — do not commit a failure you haven't
   investigated. Record any deliberately deferred failure in KNOWN_ISSUES.md
   (test name + brief reason), and clear it before declaring the branch done.
7.  A feature branch is only complete once the full test suite passes and
    KNOWN_ISSUES.md is empty for that branch.