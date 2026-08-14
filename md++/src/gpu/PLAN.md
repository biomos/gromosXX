# GROMOS CUDA Implementation — Plan

Status: **design agreed, implementation not started**. This document is the decision
record and step-by-step roadmap for the CUDA rewrite in `src/gpu`, replacing
`src/gpu_legacy`. It supersedes `notes.md` as the source of truth for *decisions*;
`notes.md` remains as historical brainstorming and is still useful for the parts
it left as open sketches (e.g. exact kernel signatures), but where it conflicts
with this document, this document wins.

Last commit before this plan was written by a human: `a5118b63d`. The commit(s)
after that (AI-generated: `cuda_nonbonded_interaction.cc`, `nb_kernels.*`,
`leap_frog_gpu.cu`, `cuda_lj_params.*`, `cuda_nb_sim_params.h` and the CPU-pairlist
+ copy-to-GPU pipeline they implement) are **to be discarded**, per §7. They are
a plausible correctness-first prototype but violate the architecture decided
here (pairlist stays on GPU, no user-facing pairlist choice under CUDA, no
per-step full-array re-upload). Salvageable pieces (LJ parameter matrix layout,
reaction-field constant precomputation, the general shape of the force kernel)
should be re-derived against the real pairlist output, not copy-pasted.

## 1. Goals and non-goals

- **Numerical correctness first.** Every GPU code path must reproduce CPU
  results within precision-appropriate tolerance before it's considered done.
  No GPU path ships without a comparison test against the CPU reference.
- **Maintainable by less-experienced contributors.** Prefer explicit, small,
  well-named classes over generic/templated machinery where both are equally
  fast. Every "clever" piece (constexpr backend dispatch, tile bitmasks,
  identity-keyed caches) needs a comment explaining *why*, and a home in this
  document.
- **Fast, not state-of-the-art.** We are not chasing GROMACS/OpenMM-level
  performance. We avoid the mistakes that make things *needlessly* slow
  (redundant host-device transfers being the big one — see §4), but we don't
  need quantization, warp-specialized kernels, or multi-GPU domain
  decomposition for v1.
- **Whole-step-on-GPU is the actual performance goal.** Individual algorithms
  are cheap; host↔device transfer is what's expensive. The target is to run
  as much of one MD step on GPU as possible so data never has to leave the
  device mid-step. This reframes several decisions below (esp. §3.1 and §8).

## 2. Decisions log

Each decision below was discussed and agreed on explicitly; treat this section
as binding unless revisited with the same rigor.

| # | Topic | Decision |
|---|-------|----------|
| D1 | GPU mirrors of `Topology`/`Configuration` | Owned by `CudaManager`, **not** by the core classes. Identity-keyed cache (§3.2). |
| D2 | `CudaManager` role | Stays the single gate to CUDA. Internals simplified: no worker-thread/queue layer, move-only types (§3.1). |
| D3 | `CudaMemoryManager` | Keep as a permanent (not debug-only) allocation-tracking layer. Fix the bugs in §5. Standardize on one smart-pointer family. |
| D4 | Precision | `FP_PRECISION` compile-time switch (1=float, 2=mixed, 3=double) stays; **default becomes mixed (2)**. |
| D5 | Pairlist algorithm choice under CUDA | **Not user-selectable.** Exactly one GPU-native tile/cell pairlist. Hard error at init if `accelerator = cuda` and the pairlist algorithm isn't the CUDA one — no automatic CPU-pairlist-then-copy fallback, ever (§6). |
| D6 | Pairlist internal design | One kernel family covering chargegroup/atomic cutoff and twin/single range, built around `TileContainerT`'s existing solute/solvent split (§6.2). |
| D7 | Verlet skin buffer | New optional `skin` field on the existing `PAIRLIST` block, default `0.0` (backward compatible). Algorithms that ignore it must warn if it's set non-zero (§6.3). |
| D8 | `make_algorithm` fallback semantics | For ordinary algorithms: CPU-only build → CPU always, error if GPU-only algorithm instantiated. GPU build + `accelerator=cuda` + algorithm has a GPU variant → **silently** use it (this is the fast, wanted path). GPU build + `accelerator=cuda` + algorithm has **no** GPU variant → fall back to CPU but **emit a warning** (this causes a hidden host-device transfer and the user should know). This warning path is not yet implemented (§8). |
| D9 | Pairlist/nonbonded vs. general algorithms | The hard-error rule (D5) is specific to pairlist/nonbonded, precisely because it's the dominant cost and a silent fallback there would look like "running on GPU" while actually eating the full CPU pairlist cost. Everywhere else, graceful fallback-with-warning (D8) is correct and must stay. |

## 3. Core infrastructure

### 3.1 CudaManager / CudaDevice

Rename `CudaDeviceManager` → `CudaDevice` (matches its role: one GPU, its
streams, its data). Drop `CudaDeviceWorker` entirely — no host worker thread,
no task queue, no condition variable. CUDA streams already provide
asynchronous execution from a single host thread; a thread-per-stream model
adds synchronization overhead and a second thread to debug for zero benefit at
GROMOS's scale (single CPU-driven `Algorithm_Sequence`, not many independent
concurrent producers).

```
CudaManager                     // facade, one per Simulation (sim.cuda())
  └─ map<device_id, CudaDevice>  // one per selected GPU
       ├─ cudaStream_t (one or two named streams: "compute", "copy")
       ├─ CudaMemoryManager
       └─ device properties
```

Make `CudaManager`, `CudaDevice` **move-only**: delete copy constructor and
copy assignment outright. The current "shallow copy, print a warning to
`std::cerr`" pattern on classes that own a `cudaStream_t` is not a style
choice, it's a live double-free / double-stream-destroy bug — two copies of
the same handle, only one of which should ever call `cudaStreamDestroy`.
Nothing in the codebase actually needs to copy these objects; `= delete` costs
nothing.

If overlap between GPU algorithms is wanted later (per notes.md's DAG/data-
dependency sketch), implement it as explicit `cudaEvent_t` dependencies
between named streams issued from the main thread — not as OS threads.

### 3.2 GPU mirrors of Topology / Configuration

`topology::Topology` and `configuration::Configuration` become GPU-agnostic
again: remove the `m_gpu` member, `get_gpu_view()`, `copy_to_gpu()` etc. from
both classes, and drop their `#include "gpu/cuda/..."` — this also removes the
circular module dependency (today gpu/ depends on topology/configuration for
the mirror-construction constructors, and topology/configuration depend on
gpu/ for the member type).

Instead, `CudaManager` owns an **identity-keyed cache**:

```cpp
class CudaManager {
  ...
  gpu::Topology::View topology_view(const topology::Topology& topo, bool force_resync = false);
  gpu::Configuration::View configuration_view(configuration::Configuration& conf, bool sync_pos_vel = true);
private:
  std::unordered_map<const topology::Topology*, gpu::Topology> m_topologies;
  std::unordered_map<const configuration::Configuration*, gpu::Configuration> m_configurations;
  const topology::Topology* m_last_topo = nullptr;       // 1-entry fast path
  gpu::Topology* m_last_topo_gpu = nullptr;
  ...
};
```

Rationale: GROMOS can have more than one `Topology`/`Configuration` object
live at once in a process (perturbation A/B setups, EDS, etc.), so a single
cached member (my first instinct) is wrong. Keying by object address costs a
hashmap lookup per call, mitigated by a 1-entry fast path for the (by far most
common) single-topology case. `Topology` data barely changes during a run
(occasional lambda/perturbation updates); `Configuration` positions/velocities
change every step — that's why `configuration_view` always takes an explicit
sync flag (mirrors `copy_pos_vel_to_device` vs. full `copy_to_device`) while
`topology_view`'s resync is opt-in and rare.

Access pattern from any algorithm: `sim.cuda().configuration_view(conf)` —
consistent with how `sim.cuda()`/`sim.mpi()`/`sim.openmp()` already expose
cross-cutting concerns through `Simulation`, which is already the fixed 3rd
argument of every `Algorithm::init/apply`. No new call-site plumbing needed.

**Cache key lifetime (must resolve before implementing, not after):** keying
purely on object address (`const topology::Topology*`) is unsafe on its own —
if a `Topology`/`Configuration` is destroyed and a new, unrelated object
happens to be allocated at the same address, the cache would silently hand
back a stale GPU mirror for the wrong object. This is worse than a crash: it
is silent wrong-data corruption, which conflicts directly with this
project's #1 goal (numerical correctness first). Removing `get_gpu_view()`/
`copy_to_gpu()` from `Topology`/`Configuration` (which is the right call, for
the header-dependency reasons above) also removes the natural place a
destructor-driven cache invalidation would have lived, so this needs an
explicit replacement, not an assumption:

- **Chosen approach:** give `Topology`/`Configuration` a plain `size_t`
  identity token member (assigned from a process-wide monotonic counter at
  construction — no GPU dependency, no new header coupling). Key the cache
  on `(pointer, token)`, or maintain a small non-GPU-aware
  `topology::Topology::id()` / `configuration::Configuration::id()` accessor
  and key on the token alone. A lookup whose stored token doesn't match the
  live object's current token is treated as a miss (rebuild the mirror),
  never as a hit against stale data.
- If, after investigation, `Topology`/`Configuration` genuinely never get
  destroyed-and-reallocated mid-run in current GROMOS usage, that's an
  acceptable reason to defer the token and rely on address-only keying — but
  that must become a written invariant with a comment at the cache
  declaration (e.g. "safe because X"), not a silent assumption, since it's
  exactly the kind of thing a later change (e.g. topology reload for REMD)
  could violate without anyone noticing.

**Thread-safety:** `CudaManager` and its caches (including the 1-entry fast
path) are **not** thread-safe and are not meant to be — all access must come
from the single CPU-driving thread, consistent with the whole-step-on-GPU /
no-worker-thread model (§3.1). Concurrent replicas or EDS legs (the case
that motivates the multi-entry map in the first place, per the rationale
above) must use separate `CudaManager` instances, never shared access to
one from multiple threads. State this as an explicit precondition/comment on
the class, not an implicit assumption.

### 3.3 CudaMemoryManager

Keep it as a real, always-on layer (not debug-only) — GPU allocations happen
at setup/resize time, not per-step, so the tracking cost is paid rarely. Its
job: device introspection, allocation bookkeeping (leak detection, "did we
actually own this pointer" checks), and a pinned-host staging-buffer pool for
async transfers (reuse buffers instead of `cudaHostAlloc`/`cudaFreeHost` every
time).

Known correctness gaps to fix before anything else builds on this class
(these are not new "issues found" so much as an acknowledged half-finished
implementation — listed here so the fix is a checklist, not a redesign):

1. `CudaMemoryManager(int device_id)` constructor never stores `m_device_id`.
2. `allocate<T>(count)` ignores `count`, always allocates `sizeof(T)`.
3. `copy_to_device(vector)` therefore overflows the allocation it just made.
4. `allocate()` returns `m_unique_ptr(devptr)` — no such member/function exists.
5. The returned `cumm_unique_ptr<T>` must be constructed with the manager
   instance in its deleter (`CuDeleter<T, CudaMemoryManager>(this)`), or
   released memory never leaves `m_allocations` and a later manager-destructor
   free can double-free.
6. `cupointer.h`'s `cu_uniquevar_ptr<T>`/`cu_uniquearr_ptr<T>` are a second,
   manager-less smart-pointer family that overlaps with `cumm_unique_ptr<T>`.
   **Decision: standardize on the manager-aware family** (`cumm_unique_ptr`);
   keep the manager-less `CuDeleter<T,void>` only for code with no manager in
   scope (should be rare/none once CudaManager is wired everywhere).
7. `cu_shared_ptr<T> = std::shared_ptr<T[]>` has no deleter — the default
   `delete[]` is undefined behavior on a `cudaMalloc`'d pointer. Either give
   it a mandatory custom deleter at every construction site or remove the
   alias until there's an actual use for shared ownership of device memory.
8. `cuda_memory_manager.tcc` is orphaned (references `allocate_device_memory`
   and a `Variable<T>` nested class that don't exist in the current header).
   Delete it.

### 3.4 Memory/allocation philosophy for algorithm authors

No central "ask the registry for a slot" API for algorithm-private GPU state.
An algorithm that needs persistent GPU scratch data just owns a
`gpu::cuvector<T>` (or a `cumm_unique_ptr<T>` from `CudaMemoryManager`) as an
ordinary class member, exactly like a CPU algorithm owns a `math::VArray`
member today. RAII ownership *is* the memory management story for per-
algorithm data; `CudaManager`/`CudaMemoryManager` exist for what's genuinely
shared (device selection, the Topology/Configuration mirrors, the pinned
staging pool) — not as a mandatory intermediary for every allocation.

### 3.5 What's already solid — keep as-is

- `gpu::cuvector<T>` / `gpu::cuhvector<T>` (`cuvector.h`) and the
  `CuMAllocator`/`CuHAllocator` allocators — correct, idiomatic
  `Allocator`-concept implementations.
- `gpu::Container<T>` (`container.h/.cu`) — the pitched 2D jagged-array
  device container, with atomic `push_back`/`reserve_strip` — solid, keep.
- `Interaction_TileT` / `TileVecT` / `TileContainerT` (`tile.h`) — the right
  shape for a GPU pairlist; not yet consumed by a real pairlist build, that's
  the work in §6.
- `gpu::Periodicity<BOUNDARY>` (`math/periodicity.h`) — including
  `prepare_chargegroup`, cell indexing via `get_cell`, and the morton-index
  hook — this is most of what §6's cell/tile build needs already.
- `HOSTDEVICE`/`DISABLED()`/`DISABLED_VOID()` macro strategy in `cuheader.h`
  — clean, keep using it for every CPU-stub method.
- The `AlgorithmB<Backend>` / `has_gpu_backend_v` / `make_algorithm<AlgT>`
  machinery in `algorithm/algorithm.h` and `util/backend.h` — validated by
  real usage in `Remove_COM_Motion`, `Lattice_Shift_Tracker`, `Leap_Frog_*`.
  Keep the pattern; see §8 for the one missing piece (the fallback warning).

### 3.6 What to delete outright

- `gpu/cuda/memory/array.h` — `CudaArray`/`CudaManagedArray`: doesn't compile
  (constructor named `CudaArray` inside `CudaManagedArray`, references
  undeclared `device_ptr`), unused, listed in `CMakeLists.txt` as a header
  only. Delete the file and its `CMakeLists.txt` entry.
- `gpu/cuda/manager/cuda_memory_manager.tcc` — see §3.3 item 8.
- `gpu/cuda/algorithm/constraints.{h,cc,cu}` — dead: `constraints.h` declares
  nothing (the real declarations are commented out), `.cc`/`.cu` define
  `gpu::algorithm::remove_com_translation/rotation` with mismatched
  signatures (`.cc` takes `const Simulation&`, `.cu` takes `Simulation&`) and
  aren't even referenced by `CMakeLists.txt` (`constraints.cu` is commented
  out of `GPU_SOURCES`). This work is superseded by the already-real
  `Remove_COM_Motion<gpuBackend>` template specialization in
  `algorithm/constraints/remove_com_motion_gpu.cc` — except that file is
  currently a byte-for-byte copy of the CPU loop with no actual GPU code (see
  §8, it's the first real algorithm to port).

### 3.7 Audit before deleting old APIs

Before removing `conf.copy_to_gpu()`, `conf.m_gpu`, `topo.get_gpu_view()`,
`conf.copy_pos_vel_to_gpu()`, or any other symbol listed for removal in §3.2/
§3.6/§7: grep the **full** source tree for each symbol name, not just the
files already enumerated in those sections. The file lists above reflect
what's known today; a mechanical migration should not rely on that list
being exhaustive. Treat any unexpected hit as a stop-and-report case rather
than silently patching around it — an unaccounted-for call site is a sign
the removal's blast radius wasn't fully understood yet.

## 4. Precision policy

Keep `precision.h`'s `FP_PRECISION` compile-time switch (1/2/3 = float/mixed/
double via `FPPolicy`), but flip the **default** in the header from
`FP_PRECISION 1` (float-only) to `FP_PRECISION 2` (mixed): bulk force/energy
math in `FPL_TYPE` (float), accumulation in `FPH_TYPE` (double). Float-only
and double-only remain available as explicit build options (debugging,
accuracy comparisons, perf comparisons) via `-DFP_PRECISION=1` or `=3`.

**Scope this explicitly before flipping it:** confirm whether `FP_PRECISION`
is consumed only by GPU-portable code (`FPL_TYPE`/`FPH_TYPE` call sites), or
whether it also changes numerics in CPU-only paths / the default CPU-only
build. This matters because §9.6 requires existing regression tests to pass
"unchanged in behavior" — that claim only holds at the bit level if the
default flip is scoped to GPU-portable code paths that don't run in a
CPU-only build. If the flip does affect CPU-only numerics too, treat it as
its own numerically-visible change: update the affected regression
baselines' tolerances explicitly and call it out in that commit's message,
rather than folding it silently into "rewrite of the accelerator layer, not
of physics." (Roadmap step 4 already lands this separately from the
pairlist/kernel work, which is the right sequencing — this note is about
making the regression-test implication explicit in the text too, not about
reordering it.)

## 5. Backend dispatch pattern (for new algorithms)

This is already implemented and validated (`util::backend.h`,
`algorithm::AlgorithmB<Backend>`, `algorithm::make_algorithm<AlgT>`). New
GPU-portable algorithms follow the existing recipe (see `notes.md`
"How to implement GPU variant of an algorithm?" for the mechanical steps —
that part of notes.md is still accurate). The one behavioral gap is the
fallback warning, tracked in §8/D8.

## 6. Pairlist and nonbonded interaction

### 6.1 The rule

Once `sim.param().gpu.accelerator == gpu_cuda`, the pairlist **must** be the
GPU-native one. This is a hard requirement, validated at init time, with no
automatic fallback:

```cpp
if (sim.param().gpu.accelerator == simulation::gpu_cuda &&
    sim.param().pairlist.grid != simulation::pl_cuda /* new selector value */) {
  io::messages.add(
    "GPU acceleration requires PAIRLIST algorithm = cuda (set PAIRLIST/GRID accordingly)",
    "CUDA_Nonbonded_Interaction", io::message::error);
  return 1;
}
```

(Illustrative: `simulation::pl_cuda` names a decision not yet made — see §11's
open question on whether the new selector is a proper `pairlist_enum` value
or reuses the existing untyped `int`. Resolve that during roadmap step 5;
this snippet shows the intended check, not a settled API.)

No "CPU pairlist + copy to GPU" path exists as a permanent feature — that's
precisely the shape of the discarded AI commit, and if it were always
available it would silently become the default under any small
misconfiguration, defeating the entire point (whole-step-on-GPU, no
transfers).

This is enforced separately from the `make_algorithm<AlgT>` mechanism,
because pairlist-algorithm selection isn't currently done through that
mechanism — it's chosen by a dedicated factory keyed on the existing
`param().pairlist.grid` integer selector among the ~6 concrete
`Pairlist_Algorithm` subclasses (`Standard_Pairlist_Algorithm`,
`..._Atomic`, `Extended_Grid_...`, `Grid_Cell_...`). Add one more concrete
value/class to that same selector (the GPU-native one) rather than bending
`make_algorithm` to a use case it wasn't shaped for.

### 6.2 One pairlist implementation, two orthogonal axes

There is exactly one GPU-native pairlist algorithm. What currently looks like
4 combinations (chargegroup×twin-range, chargegroup×single-range,
atomic×twin-range, atomic×single-range) is actually one kernel family with
two independent axes, built around the `TileContainerT` structure that
already exists (`solute_short/solute_long/solvent_short/solvent_long/
*_candidates`):

1. **Candidate build** (expensive, O(N), infrequent): build the cell/tile
   grid at `cutoff_long + skin` (Verlet buffer). Only re-run every
   `pairlist.skip_step` steps — same cadence knob that exists today. Produces
   the `*_candidates` tile lists. Port this from `gpu_legacy`'s
   `grid_pairlist.cu`/`grid.cu` (real prior art with measured perf
   characteristics — "+15% from unordered writes", "100x from the monolith
   per-cell-pair-of-blocks kernel over the looping version" — re-derive
   against that, don't re-invent from scratch), rebuilt on top of
   `gpu::Container`/`gpu::TileVecT`/`gpu::Periodicity` instead of the legacy
   `cukernel::Grid`/`Pairlist` types.
2. **Short/long classification** (cheap, O(candidates), every step): a
   distance check of each candidate tile/pair against `cutoff_short²`,
   bucketing into `solute_short`/`solute_long`/etc. When
   `cutoff_short == cutoff_long` (single-range / "fastest" mode), everything
   lands in `short` — this is the *general* case with equal cutoffs, not a
   separate code path, so single-range costs nothing extra over twin-range
   plumbing-wise.
3. **Chargegroup vs. atomic cutoff**: a template/`constexpr` parameter
   selecting which position feeds the distance test (chargegroup
   center-of-geometry via `Periodicity::prepare_chargegroup`, vs. individual
   atom position) and which unit is placed into cells. The cell/tile
   traversal is otherwise identical between the two.

Keep the existing solute/solvent split: solvent chargegroups are small and
uniform (pack into regular fixed-size tiles cheaply), solute chargegroups are
variable-size (need the more flexible candidate-list path) — this is
presumably why `TileContainerT` already separates them, and it's the right
call.

Exclusions (1-2/1-3/1-4 from the topology) get baked into the tile bitmask
(`Interaction_TileT::mask`) at candidate-build time. This is
correctness-critical and must be validated bit-exact against the CPU
pairlist's exclusion handling (§9).

### 6.3 Verlet skin parameter

New field on the existing `PAIRLIST` block (`simulation::Parameter::
plist_struct`): `double skin` (default `0.0`, backward compatible — existing
input files behave identically). Semantics: the candidate build radius is
`cutoff_long + skin`; `skin = 0` reproduces today's exact rebuild-every-
`skip_step`-with-no-buffer behavior. Any pairlist algorithm that doesn't use
`skin` (i.e. every existing CPU algorithm) must warn via `io::messages` if the
input file sets it non-zero — checked in that algorithm's own `init()`, since
it's the algorithm itself that knows whether it honors the parameter.

## 7. Discard list (the AI-generated commit)

Discard, or treat purely as reference material to be re-derived rather than
kept:

- `interaction/nonbonded/interaction/cuda_nonbonded_interaction.cc/.h` — CPU
  pairlist + full per-step position upload + one-thread-per-pair kernel.
  Violates D5/D6.
- `gpu/cuda/interaction/nonbonded/kernels/nb_kernels.{h,cu}` — the force
  kernel shape (flat `uint2` pair list, atomics for Newton's 3rd law) needs
  to become a tile-consuming kernel instead. The reaction-field math itself
  is fine and can be reused once the input format changes.
- `algorithm/integration/leap_frog_gpu.cu` — reasonable kernel bodies, but
  built on the now-removed `conf.copy_pos_vel_to_gpu()`/`conf.m_gpu` API
  (§3.2). Re-port against `sim.cuda().configuration_view(conf)`.
- `gpu/cuda/interaction/nonbonded/cuda_lj_params.{h,cu}`,
  `cuda_nb_sim_params.h` — the parameter-matrix layout and RF-constant
  precomputation are fine and can stay; they don't depend on the pairlist
  format.
- Keep `gpu/cuda/interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl_gpu.{h,cu}`
  as the file that gets *filled in* per §6 — it's currently a stub calling
  now-removed APIs (`conf.copy_to_gpu()`, `topo.get_gpu_view()`), but its
  shape (pImpl specialization of `CUDA_Pairlist_Algorithm_Impl<gpuBackend>`)
  is exactly right per D1/D2 (pImpl pattern from `notes.md`, validated at
  commit `4f78bf0`).

## 8. `make_algorithm` fallback warning (D8/D9) — the one missing behavioral piece

Current behavior (`algorithm/algorithm.h::make_algorithm`):

```cpp
template <template <typename> class AlgT, typename... Args>
Algorithm* make_algorithm(const simulation::Simulation & sim, Args&&... args) {
  if constexpr (util::has_gpu_backend_v<AlgT>) {
    if (sim.param().gpu.accelerator == simulation::gpu_cuda) {
      return new AlgT<util::gpuBackend>(std::forward<Args>(args)...);
    }
  }
  return new AlgT<util::cpuBackend>(std::forward<Args>(args)...);
}
```

This silently falls back to CPU whenever `accelerator == gpu_cuda` but
`AlgT` has no GPU backend. Per D8, that fallback should stay (it must — most
algorithms won't be ported for a long time and the run must still work), but
it should **warn**, because every silent fallback here is a hidden
host↔device round trip working against the whole-step-on-GPU goal. Add:

```cpp
template <template <typename> class AlgT, typename... Args>
Algorithm* make_algorithm(const simulation::Simulation & sim, Args&&... args) {
  if constexpr (util::has_gpu_backend_v<AlgT>) {
    if (sim.param().gpu.accelerator == simulation::gpu_cuda) {
      return new AlgT<util::gpuBackend>(std::forward<Args>(args)...);
    }
  }
  if (sim.param().gpu.accelerator == simulation::gpu_cuda) {
    io::messages.add(
      "GPU acceleration requested but no GPU backend available for this algorithm; "
      "falling back to CPU (adds a host<->device transfer)",
      "make_algorithm", io::message::warning);
  }
  return new AlgT<util::cpuBackend>(std::forward<Args>(args)...);
}
```

Apply the same to `make_unique_algorithm`. This is a small, mechanical,
low-risk first patch — good candidate for the first real commit of this
effort.

## 9. Testing and validation strategy

Numerical correctness is the top stated goal, so tests are not an
afterthought:

1. **Pairlist equivalence tests.** For a handful of representative systems
   (vacuum, rectangular PBC, triclinic; solute+solvent mix), the GPU-native
   pairlist's output (as sets of interacting atom pairs, per solute/solvent
   and short/long bucket) must exactly match `Standard_Pairlist_Algorithm`'s
   output at `skin = 0`. This includes exclusions. Run this as an actual
   automated test, not a manual spot check — it's the correctness gate
   before any GPU energy/force number can be trusted.
2. **Cross-CPU-algorithm equivalence.** While at it, add the same style of
   test between the existing CPU pairlist algorithms
   (`Standard_Pairlist_Algorithm` vs `..._Atomic` vs `Extended_Grid_...` vs
   `Grid_Cell_...`) — this doesn't exist yet and it's cheap insurance for the
   CPU side too.
3. **Force/energy comparison.** GPU nonbonded forces and energies vs. CPU,
   same input, tolerance appropriate to `FP_PRECISION` (tight for
   double-only, looser but bounded for mixed/float).
4. **Skin buffer drift test.** Verify that increasing `skin` (and therefore
   `skip_step`) doesn't change the trajectory beyond the expected
   Verlet-buffer tolerance — i.e. that no real pair is ever missed between
   rebuilds.
5. **Backend-fallback warning test.** Confirm `make_algorithm` emits the
   warning from §8 exactly when expected, and never for algorithms that do
   have a GPU backend.
6. **Build-matrix smoke tests.** CPU-only build must compile and run
   unaffected by anything under `USE_CUDA`; both must pass existing GROMOS
   regression tests unchanged in behavior (this is a rewrite of the
   accelerator layer, not of physics; see §4's note if the precision default
   flip turns out to affect CPU-only numerics).
7. **Host↔device transfer-count regression test.** Tests 1-6 verify
   correctness, but none of them verify the thing §1 states as the actual
   performance goal (minimize transfers, not per-kernel time) — a later
   change could reintroduce a hidden sync/copy without breaking any
   correctness test, silently drifting away from that goal. Instrument
   `CudaManager`/`CudaMemoryManager` with a counter incremented on every
   host↔device copy, and add a test asserting the per-MD-step count for
   fully-ported algorithms (once §10 step 11's milestone lands) doesn't
   regress upward. This is the gate for the stated goal itself, not just for
   correctness.

## 10. Implementation roadmap

Rough dependency order; each step should land as its own reviewable unit.

1. **Cleanup pass** (§3.6 discard list, §3.3 bug fixes, §8 warning patch).
   Low risk, unblocks everything else, and removes the confusing half-finished
   forks before new code goes on top of them.
2. **CudaManager simplification** (§3.1): rename, drop worker/queue, move-only.
3. **GPU mirror relocation** (§3.2): move `gpu::Topology`/`gpu::Configuration`
   out of the core classes into `CudaManager`'s identity-keyed cache; update
   the one real consumer so far (`Leap_Frog_*<gpuBackend>`, once re-ported)
   to the new access pattern.
4. **Precision default flip** (§4): one-line change + rebuild/test matrix.
5. **Port the legacy grid/cell pairlist** onto `Container`/`TileVecT`/
   `Periodicity` (§6.2, step 1: candidate build). This is the largest single
   piece of new work.
6. **Short/long classification + chargegroup/atomic axis** (§6.2, steps 2-3).
7. **Pairlist equivalence tests** (§9.1-9.2) — before writing the force
   kernel, so the pairlist is trusted on its own.
8. **Force kernel rewritten against tiles** (adapting the RF math from the
   discarded `nb_kernels.cu`, consuming `Interaction_Tile` instead of a flat
   pair array).
9. **`CUDA_Nonbonded_Interaction` wired end-to-end**, hard-error gate (§6.1)
   in place.
10. **Force/energy comparison tests** (§9.3), skin drift test (§9.4).
11. **DONE.** **Re-port `Leap_Frog_*<gpuBackend>`** against `sim.cuda().
    configuration_view()` instead of the removed `conf.m_gpu` API — first
    real "runs entirely on GPU across two algorithms without a round trip"
    milestone. Landed in two parts: the re-port itself (new
    `gpu/cuda/algorithm/integration/leap_frog_kernels.{h,cu}` +
    `algorithm/integration/leap_frog_gpu.cc`, validated by
    `leap_frog_gpu.t.cc`), then a correctness fix for a real bug the re-port
    shipped with — `Leap_Frog_Position<gpuBackend>` silently used a stale
    GPU-resident velocity whenever CPU-only temperature coupling
    (`multibath.couple`) ran between it and `Leap_Frog_Velocity<gpuBackend>`.
    The fix generalizes past this one case: `CudaManager` now tracks
    per-field freshness on the `Configuration` mirror
    (`gpu_fresh_fields`/`gpu_dirty_fields`, `gpu/mirror_fields.h`) instead of
    relying on hand-picked `sync_pos_vel`/`full_resync` booleans at each call
    site, with `Algorithm_Sequence::run()` centrally flushing/invalidating
    it around every algorithm's `apply()` via a new
    `Algorithm::gpu_mirror_touches()` hook (default-safe for every existing
    CPU algorithm; see `TILE_PAIRLIST_DESIGN.md`'s "Follow-up: GPU-mirror
    freshness tracking" section for the full design and the regression test
    that caught it).
12. **Shake/constraints on GPU** (currently only stubs in
    `gpu/cuda/algorithm/constraints.*`, which are being deleted per §3.6 —
    this is new work using the same backend-dispatch pattern as
    `Remove_COM_Motion`).
13. Broaden algorithm coverage incrementally, each addition removing one more
    host↔device round trip per step, each with a fallback warning (§8) until
    ported. **In progress:** `Remove_COM_Motion<gpuBackend>` -- done. It
    already had the `Backend` template slot and was already selected by
    `make_algorithm` under `accelerator = cuda`, but its `apply()` was the
    plain CPU loop wearing a `gpuBackend` label (no kernel at all) --
    replaced with real reductions (translation: sum(m·v); rotation:
    angular momentum + inertia tensor, same two-pass shape) in
    `gpu/cuda/algorithm/constraints/remove_com_motion_kernels.{h,cu}`,
    matching `remove_com_motion_cpu.cc`'s formulas exactly. Verified against
    the CPU reference (`remove_com_motion_gpu.t.cc`) and with zero
    `compute-sanitizer --tool memcheck` errors.

    **Also done:** `Temperature_Calculation<gpuBackend>` +
    `Berendsen_Thermostat<gpuBackend>`. Unlike `Remove_COM_Motion`, neither
    class had a `Backend` template at all before this -- both converted
    from plain classes to templates (`create_md_sequence.cc` already
    called both through `make_algorithm<T>`, so no call-site changes were
    needed once they became templates; it was silently resolving to
    `make_algorithm`'s legacy non-templated overload before). The real
    complexity: both need a temperature group's centre-of-mass velocity,
    an inherently two-level (atom → molecule) reduction, done via a new
    `launch_group_velocity_reduce()` kernel
    (`gpu/cuda/algorithm/temperature/temperature_kernels.{h,cu}`) bucketed
    by a per-atom temperature-group index (built by hand from
    `topo.temperature_groups()`, since -- found the hard way -- that CSR
    list uses a *different* convention than `energy_groups()`: `num_groups
    + 1` entries with a leading `0` and exclusive-end boundaries, not
    `num_groups` entries with an inclusive end and an implicit start at 0;
    confirmed directly from `in_topology.cc`'s parsing code, not
    guessed). `Berendsen_Thermostat<gpuBackend>` leaves the scaled
    velocity resident on the GPU mirror (`mark_gpu_dirty()`, no
    sync-back) -- it sits between `Leap_Frog_Velocity<gpuBackend>` and
    `Leap_Frog_Position<gpuBackend>` in `create_md_sequence.cc`, so that
    whole three-algorithm chain now runs with zero host round trips when
    temperature coupling is on, the exact scenario this session's
    GPU-mirror freshness-tracking fix was built around. Verified against
    the CPU reference on a synthetic 3-bath setup exercising both of
    `Thermostat::scale()`'s cases (`com_bath == ir_bath` and `!=`,
    `temperature_gpu.t.cc`) and with zero `compute-sanitizer` errors.

    **Also done:** `CUDA_Quartic_Bond_Interaction` (bonded forces, first
    of four planned terms). Unlike the `Algorithm`/`Backend` template
    pattern used above, bonded terms use the `interaction::Interaction`
    dispatch mechanism (same family as `CUDA_Nonbonded_Interaction`):
    `create_bonded.cc` now instantiates `CUDA_Quartic_Bond_Interaction`
    instead of `Quartic_Bond_Interaction` when `accelerator == gpu_cuda`
    and perturbation is off (perturbation is hard-errored in the new
    class's own `init()`, same convention as every other CUDA gate).
    Structurally much simpler than nonbonded: no pairlist at all -- the
    term list (atom-index pairs + type index, from
    `topo.solute().bonds()`/`topo.bond_types_quart()`) is static and
    small (aladip: ~80 bonds), uploaded once in `init()`, consumed by a
    plain one-thread-per-term kernel
    (`gpu/cuda/interaction/bonded/quartic_bond_kernels.{h,cu}`) with
    direct global `atomicAdd`s (no shared-memory bucketing -- term counts
    are too small to bother). Virial is unconditional, matching the CPU
    source's dead `if (V == math::atomic_virial)` gate. **Precision
    finding, worth remembering for every future bonded term:** the
    default `FP_PRECISION=1` (float-only) build stores the position
    mirror in `FPL_TYPE` (float), and `dist2 - r0^2` is a near-
    cancellation for a bond near its equilibrium length -- aladip's
    largest bond force constant (~1.57e7) turns the position mirror's
    ordinary float-truncation error into a force discrepancy an order of
    magnitude above the `1e-4` relative tolerance used for nonbonded.
    Confirmed by direct calculation (matches observed magnitude) and by
    testing that computing the subtraction itself in `double` inside the
    kernel made no difference (the error is baked into the input
    position's float representation, not the arithmetic) -- this is an
    inherent property of the float position mirror, not a kernel bug, so
    every future bonded term's test should budget for it (`5e-3` relative
    used in `quartic_bond_gpu.t.cc`) rather than reusing nonbonded's
    tighter tolerance verbatim. Verified against the CPU reference
    (`quartic_bond_gpu.t.cc`) and zero `compute-sanitizer` errors.

    Still CPU-only: Angle/Improper Dihedral/Torsional Dihedral bonded
    terms, constraints (SHAKE/SETTLE/LINCS -- step 12 above),
    `NoseHoover_Thermostat`, `Berendsen_Barostat`, `Pressure_Calculation`.

## 11. Open questions (revisit later, not blocking the plan above)

- Multi-GPU load balancing / domain decomposition — out of scope for v1
  (notes.md's sketch stands as a future direction, not a current commitment).
- Whether the pinned-host staging pool (§3.3) needs to be sized/tunable via
  input parameters, or a fixed internal default is fine to start.
- Where exactly the new `pairlist.grid` selector value/constant should live
  (a new `simulation::pairlist_enum` vs. reusing the existing untyped `int` —
  worth cleaning up either way, decide during step 5 of the roadmap).
- Multiple concurrent `CudaManager` instances (e.g. replica exchange, EDS
  legs) potentially sharing one physical GPU: §3.2's per-`Simulation`
  `CudaManager` model assumes single-threaded access per manager (see
  §3.2's thread-safety note), but doesn't yet say how multiple managers
  sharing one device should coexist — separate CUDA contexts per process is
  the likely answer if replicas are separate processes, but this hasn't
  been confirmed against how GROMOS actually runs REMD/EDS today. Revisit
  before any multi-replica GPU deployment.
