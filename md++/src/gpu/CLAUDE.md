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

- Always work on a branch named `claude/<short-task-description>`, never commit directly to `main`
- Commit after each logical, working change with a clear message
- Never force-push
- Never merge to `main` yourself — leave that for human review
- Run tests before every commit; do not commit code that fails tests

## Build

Build out-of-tree in a dedicated directory so source stays clean:

​```bash
mkdir -p build && cd build
../Config.sh          # only needed once, or after configure.ac changes
../configure --disable-debug --disable-shared --enable-static --enable-openmp --with-cuda
make -j $(nproc)
make check
​```

Rules:
- Do not run `Config.sh` on every build — only if `configure` is missing or `configure.ac`/`Makefile.am` changed.
- Always run from a clean `build/` directory if you changed `configure` flags; run `make clean` first if switching flags on an existing build dir.
- A build is only "done" when `make -j $(nproc)` exits 0 AND `make check` exits 0. Do not report success otherwise.

## Build/test failure procedure

If `configure` fails:
1. Read the exact error in the output (not just the summary at the end).
2. Check `config.log` for the specific failing check.
3. Do not disable a feature to work around a missing dependency — report what's missing instead, unless explicitly told the dependency is optional.

If `make` fails:
1. Fix the actual compile/link error shown — do not suppress warnings-as-errors by weakening flags unless the warning is a false positive you can justify in the commit message.
2. Re-run `make -j $(nproc)` after each fix; do not batch multiple unverified fixes before re-testing.

If `make check` fails:
1. Identify which specific test(s) failed and read their output/log (usually in `build/*.log` or `test-suite.log`).
2. Distinguish a real bug (code produces wrong result) from a broken test (test itself is outdated/incorrect) — fix the actual bug by default; only modify a test if you can justify the test itself is wrong.
3. Re-run only the failed test first if possible, then the full suite, before considering it fixed.
4. If a fix is non-obvious after 2-3 attempts, stop and leave a clear note in the commit/PR description rather than making increasingly speculative changes.

## Writing new tests

- New tests go in src/check
- Follow the existing test file naming/structure — look at 2-3 existing tests before writing a new one
- Register new test files in the relevant `Makefile.am` (`TESTS = ...` list) so `make check` picks them up automatically
- Every bug fix should come with a regression test that fails before the fix and passes after, where practical
- Run `make check` after adding a test to confirm it's actually being picked up by the harness (not silently skipped)

## CUDA build notes
- Verify `nvcc --version` succeeds before building; if CUDA isn't found, do not silently drop `--with-cuda` — report it.
- If a specific GPU compute capability / arch flag is required and not autodetected, check `configure --help` output for a `--with-cuda-arch`-style flag rather than guessing.

## Standard workflow
1. Make code change
2. Build (`make -j $(nproc)`)
3. Test (`make check`)
4. If failures: diagnose → fix → rebuild → retest (do not skip straight to committing)
5. Commit only once build + full test suite pass