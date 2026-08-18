# Real tile-based GPU pairlist: design

Status: **design proposal, not yet implemented.** This is `PLAN.md` §10
roadmap step 5 ("port the legacy grid/cell pairlist onto
`Container`/`TileVecT`/`Periodicity`" — "the largest single piece of new
work"). `PLAN.md` D5/D6 and §6.2 already fixed the high-level shape (one
GPU-native pairlist, candidate-build + short/long-classify as two phases,
chargegroup/atomic as a template axis); this document works out the
concrete kernels, data structures, and gaps needed to actually build it,
and flags the sub-decisions that still need sign-off before code gets
written.

## 1. Inventory: what already exists and what's missing

**Usable as-is:**
- `gpu::Interaction_TileT<32,32>` / `TileVecT` / `TileContainerT` (tile.h) —
  the tile type and the atomic-`push_back`, overflow-tracked unified-memory
  vector holding them. `TileVecT::push_back` is finished and correct.
- `gpu::Container<T>` (container.h) — pitched jagged array with atomic
  `push_back`/`reserve_strip`. Good fit for a fixed-capacity cell list.
- `gpu::Periodicity<BOUNDARY>::set_cell_size`/`get_cell`/
  `prepare_chargegroup` (math/periodicity.h) — cell sizing and per-
  chargegroup cell assignment + cog computation + box-wrapping.
  **Gap: `set_cell_size` only has real logic for `vacuum` and
  `rectangular`; `triclinic`/`truncoct` fall into an unimplemented branch
  (prints to `std::cerr`, leaves cell size zero).** See §4.
- `gpu::prepare_cog_kernel`/`put_chargegroups_into_box_kernel`
  (kernels/periodicity.cu) — real, mostly-working: computes cog, wraps into
  box, computes cell id (`ushort4` with a Morton index in `.w`) per
  chargegroup, explicitly instantiated for vacuum/rectangular/triclinic
  (triclinic instantiated but not functionally correct per the gap above).
  **Bug to fix while touching this: `CUDA_Pairlist_Algorithm_Impl_gpu.cu`'s
  `_prepare_cog` launches it with `dim3 dimGrid(1)` — a single block
  regardless of chargegroup count.**

**Real prior art, not reusable verbatim (`gpu_legacy`):**
- `cukernel::interaction::GridT<T>` (grid.h/grid.cu) — 3-pass cell-list
  build (count via atomic increment, then an **O(cells) serial scan** per
  thread to compute each cell's start pointer, then scatter). The middle
  pass is the wrong instinct for a from-scratch implementation --
  `gpu::Container<T>`'s atomic `reserve_strip`/`push_back` gets the same
  result in one pass without a prefix-sum-shaped placeholder. Also
  contains a literal syntax error (`cuinteraction::::assign_cells`,
  double-colon) confirming this file doesn't currently compile — read for
  the algorithm shape, not the code.
- `cukernel::interaction::find_pairs_neighbour`/`find_pairs_self`
  (grid_pairlist.h/.cu) — real, documented (measured "+15% unordered",
  "100x monolith-over-looping") cell-pair distance-test kernels. Useful
  reference for the shared-memory staging pattern and the atomic
  `reserve_strip` offset-reservation idea, **but it produces a flat,
  per-atom pairlist (conceptually the same ragged shape as CPU `Pairlist`,
  just built on GPU), not tiles, and has no exclusion handling at all.**
  It answers "how do you find neighbor pairs from a cell grid on GPU," not
  "how do you fill an `Interaction_TileT`."

**Missing entirely (new work, not stubs to fill in):**
1. **Exclusions are not represented on the GPU side at all.**
   `gpu::TopologyView` today has `iac`/`mass`/`inverse_mass`/`charge`/
   `chargegroup` and counts — nothing for 1-2/1-3/1-4 exclusions.
   `topology::Topology::exclusion(i)` returns a sorted `std::vector<int>`
   of `j > i` partners, host-side, ragged, built once at topology-read time
   and static for the rest of a normal (non-perturbed-topology-reload) run.
   This needs a GPU-friendly precomputed form, built once (not per pairlist
   rebuild), that a tile-classification kernel can query to bake exclusion
   bits into `Interaction_TileT::mask`.
2. **The `skin` parameter (`PLAN.md` D7) doesn't exist on `plist_struct`
   yet.** Needed as the candidate-build radius offset
   (`cutoff_long + skin`).
3. **No atom-level (as opposed to chargegroup-level) cell assignment
   exists**, needed for the `atomic_cutoff` axis of D6.
4. **No cluster/block formation.** See §2 -- this is the actual "how do
   32 atoms become a tile row" question, and nothing in the tree answers
   it yet.

## 2. The central design decision: how do 32 atoms become a tile row?

`Interaction_TileT<32,32>` needs its 32-atom "row" and 32-atom "column" to
be *fixed-size, contiguous* blocks of *some* ordering of atoms/chargegroups
-- a cell from `Periodicity::get_cell` is not fixed-size (a cell can have 3
atoms or 80, depending on local density). Two ways to reconcile that:

- **(A) Cell-native, padded tiles.** Keep cells as the unit; when a cell
  has fewer than 32 members, pad the tile with sentinel/invalid entries
  (masked off). Tiles get generated per cell-pair (self + neighbor cells,
  as in `gpu_legacy`'s `find_pairs_self`/`find_pairs_neighbour`). Simple
  mentally, closely follows the legacy code's structure, but wastes
  occupancy whenever cells are smaller than 32 (common unless cell size is
  tuned so each cell holds close to 32 atoms, which fights against
  "cell size just above the cutoff" from `set_cell_size`), and produces a
  variable, data-dependent number of tiles per cell-pair that's awkward to
  reason about for the short/long classification pass.
- **(B) Sorted fixed-size blocks ("cluster list"), GROMACS/OpenMM-style.**
  Reorder chargegroups (or atoms, for `atomic_cutoff`) once per rebuild
  into a single spatially-coherent sequence (sort by cell id, using the
  Morton index already computed in `ushort4.w` as the sort key -- so
  spatially close atoms end up index-adjacent), then simply chunk that
  sequence into fixed 32-wide blocks. A block's bounding sphere/box is
  compared against neighboring blocks' bounding volumes (not full N×N)
  to decide which block-pairs become candidate tiles. This is the
  standard, well-proven approach (GROMACS "nbnxm" cluster pairlist, OpenMM
  similarly) specifically because it fits fixed-width SIMD/warp tiles
  without padding waste, and it makes "how many tiles total" driven by
  block-bounding-volume tests rather than raw cell occupancy.

**Recommendation: (B).** It's the reason `Interaction_TileT` is fixed at
32×32 in the first place (matches one warp) -- committing to that shape
and then padding cell-native tiles fights the format instead of using it.
It costs one extra pass (the sort/reorder), which `gpu::Periodicity`
already computes the sort key for (Morton index), and it directly reuses
`prepare_chargegroup`'s existing cog+cell computation.

## 3. Proposed pipeline (candidate build, D6 phase 1)

**Correction (found while implementing step 5, confirmed against
`nonbonded_set.cc`):** the entire pipeline below -- candidate build *and*
short/long classification -- runs only at pairlist rebuild time (every
`pairlist.skip_step` steps), never every step. This is not an optimization
choice, it's what GROMOS's twin-range scheme actually does: which pairs are
short vs. long is decided once per rebuild; only the *force values*
recomputed from `solute_short`/`solvent_short` happen every step
(`m_outerloop.lj_crf_outerloop(..., false, ...)`, called unconditionally),
while `solute_long`/`solvent_long` forces are computed once per rebuild
(inside the `if (pairlist_update)` block) and implicitly held as a frozen
background contribution until the next rebuild. An earlier draft of this
section said classification runs "every step" -- that would build a
different (not equivalent) scheme from what `PLAN.md` §9.1's equivalence
test needs to match. Both `CUDA_Pairlist_Algorithm_Impl::reorder()` and
`build_candidates()`/`classify_tiles()` all run inside
`CUDA_Pairlist_Algorithm::update()` (the `skip_step`-gated call), not
`prepare()` (called every step, box-wrapping/cog only).

**Second correction (also found while implementing step 5):** blocks must
be **atom-indexed**, not chargegroup-indexed. `Interaction_TileT::mask` is
a 32×32 *atom*-pair bitmask; exclusions are defined between specific atom
pairs, not chargegroup pairs, so a chargegroup-indexed tile has no
granularity to mask exclusions against. Chargegroup-cutoff mode still
decides short/long at the chargegroup-COG level (matching
`Standard_Pairlist_Algorithm` exactly) -- but that's a per-atom-pair
*lookup* the classification kernel does (via each atom's owning
chargegroup), not a property of how blocks are built. So: blocks are 32
atoms each, atom bounding "position" is just that atom's own (already
box-wrapped) position, and the classification kernel looks up each atom's
owning chargegroup (binary search over `TopologyView::chargegroup`, the
existing offset array) to fetch that chargegroup's cog (`m_cg_cog`,
already populated for every chargegroup, solute and solvent, since
`Periodicity::prepare_chargegroup`'s guard was relaxed in step 3) for the
actual cutoff test. This unifies chargegroup/atomic cutoff exactly as D6
originally intended: same tiles, same kernel, only the distance-test input
(chargegroup cog-cog vs. atom-atom) changes.

Per rebuild:

1. **Cell assignment** (`prepare_cog_kernel`, already exists): computes
   each chargegroup's cog, box-wraps its atoms, and emits a Morton sort
   key. Every atom in a chargegroup inherits that chargegroup's sort key
   (not its own recomputed one) -- keeps atoms of the same chargegroup
   spatially adjacent after sorting, which is what actually matters for
   bounding-sphere tightness; a per-atom key would be no more correct,
   just more code.
2. **Sort by cell/Morton key** into a permutation array, at atom
   granularity: `(atom_sort_key[a], a)` pairs sorted via
   `thrust::sort_by_key`, separately for the solute atom range
   `[0, num_solute_atoms)` and solvent atom range
   `[num_solute_atoms, num_atoms)`.
3. **Block formation**: chunk each sorted atom order into fixed 32-wide
   blocks (`num_blocks = ceil(num_atoms_in_group / 32)`; the last block is
   simply shorter, no padding sentinel needed since every kernel computes
   the valid range directly). Compute each block's bounding sphere (center
   + radius) directly from the 32 atoms' own positions.
4. **Block-pair candidate search**: `O(num_blocks_a * num_blocks_b)` per
   call, no neighbor-cell pruning yet (deliberate v1 simplification, see
   the implementation's own comment for why this is correctness-preserving
   even though it costs more comparisons than necessary). Test bounding-
   sphere distance against `cutoff_long + skin`; survivors get pushed into
   `solute_candidates`/`solvent_candidates` as an `Interaction_TileT` with
   `index` encoding the two block ids and an all-zero mask (mask gets
   filled in step 5, not here).
5. **Exclusion + short/long classification**: one CUDA block per candidate
   tile (32×32 threads, one warp per tile row). Each thread resolves its
   atom pair, skips padding/self/duplicate-diagonal/same-chargegroup
   pairs, looks up the exclusion CSR (skip if excluded), then the owning
   chargegroups' cog-cog distance decides short/long (or out-of-range
   entirely, since the bounding-sphere test in step 4 is conservative).
   Each warp uses `__ballot_sync` to build its row's 32-bit mask word with
   no atomics; the block then pushes a new tile (same index, computed
   mask) into `solute_short`/`solute_long`/`solvent_short`/`solvent_long`
   if that mask is non-zero. `cutoff_short == cutoff_long` (single-range)
   is the degenerate case where every non-excluded, in-range pair lands in
   `*_short` and `*_long` never gets a hit, exactly as `PLAN.md` §6.2
   specifies -- no separate code path.

## 4. Scope decisions needed before implementation (not yet decided)

1. **Boundary condition coverage for v1.** `Periodicity<BOUNDARY>::
   set_cell_size` only has real logic for `vacuum`/`rectangular`.
   Triclinic cell lists need a skewed-cell nearest-image scheme (GROMACS's
   own triclinic Verlet-list code is non-trivial for exactly this reason).
   Recommendation: **v1 supports vacuum + rectangular only**, and the
   pairlist's `init()` hard-errors (`io::messages`, `io::message::error`)
   if `accelerator == gpu_cuda` and the boundary is triclinic/truncoct,
   rather than silently producing wrong cell assignments. Real triclinic
   support becomes a follow-up roadmap item once vacuum/rectangular is
   validated end-to-end.
2. **Chargegroup-cutoff first, or both axes at once?** D6 wants both
   eventually. Recommendation: implement chargegroup-cutoff first
   (`prepare_chargegroup` already exists and is closer to done), get it
   passing the equivalence test against `Standard_Pairlist_Algorithm`
   (`PLAN.md` §9.1), *then* add the atomic-cutoff template branch as a
   second, smaller, separately-reviewable step -- not both in one pass.
3. **`skin` parameter addition.** Add `double skin = 0.0` to
   `plist_struct` now (small, mechanical, per D7) so the candidate-build
   radius is `cutoff_long + skin` from the start, even though no
   displacement-tracking/drift logic (`PLAN.md` §9.4's skin-buffer-drift
   test) will exist yet -- rebuild cadence stays governed by the existing
   `skip_step` counter only. Flagging so it's clear "skin support" here
   means "the candidate radius has a buffer," not "we track whether the
   buffer was actually enough between rebuilds" -- that's real future work
   this doesn't claim to do.
4. **Exclusion upload timing.** Build the GPU exclusion structure once
   in `CUDA_Pairlist_Algorithm::init()` (topology is static for a normal
   run), not on every `update()`. If perturbation ever needs a different
   exclusion set mid-run, that's out of scope here and should hard-error
   for now rather than silently using stale exclusions.

## 5. Testing (ties to `PLAN.md` §9.1)

Before any force/energy number from this pairlist is trusted:
- Build a small standalone test (in `src/check`, per `CLAUDE.md`'s testing
  conventions) that runs both `CUDA_Pairlist_Algorithm` and
  `Standard_Pairlist_Algorithm` on the same topology/configuration
  (vacuum and rectangular, solute+solvent mix, at `skin = 0`) and asserts
  the resulting pair sets are identical -- per atom, per short/long
  bucket, including exclusions. This is the gate mentioned in `PLAN.md`
  §9.1; it doesn't exist yet for either the CPU-vs-CPU or CPU-vs-GPU case.
- Only after that test passes does force/energy comparison (`PLAN.md`
  §9.3) become meaningful to run.

## 6. Proposed implementation sequence (separately reviewable steps)

1. **Done.** Add `skin` to `plist_struct` (+ parser, + warning-if-nonzero-
   and-ignored in the four existing CPU algorithms per D7).
2. **Done.** Build the GPU exclusion structure (gap §1.1), wired into
   `gpu::TopologyView`/`gpu::Topology`, built once in the constructor.
3. **Done.** Fix `prepare_cog_kernel`'s launch config bug; extend it to
   also emit the sort key it already computes.
4. **Done.** Cell/block build: sort-by-key (Thrust) + block-bounding-
   sphere kernel + block-pair candidate kernel, chargegroup-cutoff only,
   vacuum/rectangular only, writing into `solute_candidates`/
   `solvent_candidates`. Reworked to atom-indexed blocks partway through
   (see this document's "second correction" above).
5. **Done.** Short/long classification + exclusion-mask kernel,
   chargegroup-cutoff only, writing into `solute_short`/`solute_long`/
   `solvent_short`/`solvent_long`.
6. **Done.** The pairlist-equivalence test (§5) against
   `Standard_Pairlist_Algorithm`: `src/check/pairlist_cuda_equivalence.t.cc`,
   vacuum + rectangular, both passing exactly (144/144 atoms, all four
   buckets). Found and fixed two real bugs along the way, not just design
   gaps -- both worth remembering since they're the kind of thing that
   silently gives plausible-looking wrong answers without a test like
   this one:
   - `TileVecT` owns its memory (destructor calls `cudaFree`) but was
     being passed **by value** into every kernel. Kernel launch syntax
     constructs a host-side temporary for a by-value argument, and that
     temporary's destructor runs synchronously right after the
     (asynchronous) launch is enqueued -- freeing the shared buffers
     while the kernel may still be running, and leaving the "original"
     object (e.g. `m_tiles.solute_candidates`) with dangling pointers.
     Fixed by adding `TileVecT::View` (non-owning, no destructor,
     mirrors the existing `cuvector`/`CuVArray` View pattern) and
     changing every kernel signature to take `View` instead of the
     owning `TileVecT`.
   - `estimate_candidate_capacity`'s `static_cast<unsigned>(per_block *
     num_blocks_a)` truncated instead of rounding up, silently
     underestimating capacity by up to 1 whenever the fractional part
     was large (e.g. 0.93 truncates to 0). The resulting overflow was
     real and correctly detected by `TileVecT::was_overflowed()`/
     `check_candidate_overflow`'s `io::messages.add(..., error)` call --
     but `io::messages.add` only queues a message, and the equivalence
     test wasn't calling `io::messages.display()` after running the
     algorithms, so the error was queued and never shown. Fixed both:
     `std::ceil()` instead of a truncating cast, and the test now flushes
     and checks `io::messages` after every run (`flush_messages()`),
     treating any error-or-worse severity as a hard failure rather than
     silently continuing.
7. **Done.** Atomic-cutoff axis, same kernels/tiles: `classify_tiles_kernel`
   gained a compile-time `ATOMIC_CUTOFF` template bool (alongside its
   existing `BOUNDARY` template) that switches the distance-test input
   (atom-atom `nearest_image` vs. chargegroup cog-cog) and the
   same-chargegroup handling (real distance test in atomic mode, matching
   `Standard_Pairlist_Algorithm::update_atomic`'s "no special-casing"; an
   assumed-always-short shortcut in chargegroup mode, matching
   `_update_cg`'s direct push to `solute_short`) -- both gated by an
   exclusion check first, and both skipping solvent same-chargegroup
   (same-molecule) pairs unconditionally, since those have no CSR entries
   to check either way. `reorder()` gained a new `atomic_cutoff`-only path
   (`gpu::atom_cell_kernel`, `_atom_sort_key_atomic`) that computes each
   atom's own Morton sort key from a locally box-wrapped position copy,
   since `prepare_cog()` skips the chargegroup cog/cell build entirely in
   this mode (matching `update_atomic`'s CPU reference, which never wraps
   positions either -- `nearest_image` alone handles PBC). `prepare_cog()`
   now always refreshes the GPU position mirror (`conf.copy_to_gpu()`)
   regardless of mode, since atomic-cutoff has no other call site that
   does. `CUDA_Pairlist_Algorithm::init()`'s former atomic_cutoff
   hard-error, and `prepare()`/`update()`'s runtime re-check guards against
   a toggle mid-run, are both removed -- no longer needed now that every
   stage re-reads `sim.param().pairlist.atomic_cutoff` fresh and both
   modes are real, safe-to-toggle-between implementations.

   Found and fixed a real, pre-existing correctness bug in step 5's
   chargegroup-cutoff code while building this, not just new work: the
   original `classify_tiles_kernel` unconditionally skipped *every*
   same-chargegroup pair (`if (cg1 != cg2)`), but
   `Standard_Pairlist_Algorithm::_update_cg` only does that for *solvent*
   chargegroups (structural, since solvent atoms have no exclusion-CSR
   entries) -- solute same-chargegroup pairs are supposed to go through
   a real exclusion check and, if not excluded, land in `solute_short`
   unconditionally (no distance test). The blanket skip happened to match
   `aladip`'s topology (apparently every intra-chargegroup solute pair is
   already 1-2/1-3-excluded there) so step 6's equivalence test didn't
   catch it, but it would silently drop real short-range interactions for
   any topology with a chargegroup containing a non-excluded internal
   pair. Fixed as part of this step's kernel rework (see above); the
   existing vacuum/rectangular chargegroup-cutoff equivalence cases still
   pass exactly after the fix, confirming it didn't change `aladip`'s
   result.

   Second equivalence test: extended `pairlist_cuda_equivalence.t.cc`
   (not a separate file -- `Standard_Pairlist_Algorithm` already dispatches
   chargegroup/atomic internally via the same `atomic_cutoff` flag, no
   separate "`_Atomic`" class exists) with two more cases (vacuum +
   rectangular, `atomic_cutoff = true`), 4 cases total, all passing
   exactly.

   **Follow-up fix (much later, KNOWN_ISSUES.md's "CUDA context
   corruption after runtime atomic_cutoff toggle" entry):**
   `prepare_cog()` sized `m_cg_cog` to `num_solute_chargegroups()`, but
   `prepare_cog_kernel` writes `cg_cog[cg_i]` for every chargegroup
   (solute and solvent) and `classify_tiles()` reads it back for
   solvent candidates too -- an out-of-bounds write on every
   `prepare_cog()` call whenever any solvent chargegroups exist (i.e.
   essentially always), found via `compute-sanitizer --tool memcheck`
   rather than any visible test failure (undefined behaviour that
   happened not to corrupt anything load-bearing on most runs). Fixed
   by sizing `m_cg_cog` to `num_chargegroups()`. New regression test,
   `cuda_atomic_cutoff_toggle.t.cc` (a non-perturbed system, since
   `aladip_cuda.in`'s `PERTURBATION` setting prevents `ctest` from ever
   reaching this code at all) -- verified with zero errors under
   `compute-sanitizer` across 10 toggle cycles. Worth periodically
   re-running the CUDA suite under `compute-sanitizer` even absent a
   visible failure; this bug had been silently present since step 3-5's
   original chargegroup-cutoff implementation.

8. **Done (kernel only; wiring is step 9, not started).** `PLAN.md` §10
   step 8: the LJ + reaction-field force/energy kernel, rewritten against
   tiles instead of the discarded flat-pair-array kernel
   (`gpu/cuda/interaction/nonbonded/kernels/nb_kernels.{h,cu}`, deleted --
   never actually wired into anything, per `CLAUDE.md`'s "salvageable math,
   not the surrounding data flow"). New home:
   `gpu/cuda/interaction/nonbonded/kernels/lj_crf_tiles.{h,cu}`, plus the
   two small reused-as-is pieces that were already sitting unregistered in
   the tree (`cuda_lj_params.{h,cu}`: GPU LJ parameter matrix, built once
   from `interaction::Nonbonded_Parameter`; `cuda_nb_sim_params.h`: the
   `four_pi_eps_i`/`crf_2cut3i`/`crf_cut`/cutoff constants struct) -- both
   now actually registered in `gpu/CMakeLists.txt` and used.

   `lj_crf_tile_kernel<BOUNDARY>` mirrors `classify_tiles_kernel`'s layout
   exactly (one CUDA block per tile, 32x32 threads, one warp per row,
   `row_order`/`col_other_order` convention reused verbatim) and computes
   `interaction::Nonbonded_Term::lj_crf_interaction`'s default case
   (`eps=0`, `coulomb_scaling=1`) per active mask bit: forces via
   `atomicAdd` (one thread per pair), energies via a warp-shuffle-then-
   shared-memory reduction to one `atomicAdd` per tile (into a `double`
   accumulator regardless of `FPL_TYPE`, so many small per-tile
   contributions don't lose precision the way thousands of individual
   float atomicAdds would).

   `lj_crf_tiles.h` deliberately does **not** declare the `__global__`
   kernel template or include `gpu/cuda/math/periodicity.h`: that header
   isn't host-compilable (`Periodicity::set_cell_size`/`get_cell` use bare
   `min`/`max`, which only resolve under `nvcc`), so a header a plain
   `.cc` test needs to include can't pull it in. Only a host-callable
   `launch_lj_crf_tiles(..., math::boundary_enum boundary, math::Box box,
   ...)` is public; it dispatches to the right `Periodicity<BOUNDARY>`
   instantiation internally, entirely inside the `.cu` file. Found the
   same problem transitively through `block_pairlist.h` (declares
   `__global__` kernels taking `Periodicity<BOUNDARY>` too) while writing
   the test below -- neither header is meant to be included from a `.cc`
   translation unit, only from other `.cu` files.

   Also added `TileVecT::resize(unsigned)` (`tile.h`): both `push_back`
   overloads are `__device__`-only (atomic-append from a kernel), so
   there was no way to populate a `TileVecT`'s `size()` from host code at
   all -- needed for the test below to hand-build one tile. Calling the
   `__device__`-only `push_back` from host code doesn't fail to compile
   under a plain (non-`nvcc`) `.cc` compiler (attributes just vanish), it
   silently links to nothing usable and crashes at runtime with no
   diagnostic -- found and fixed while building the test, not a
   pre-existing bug (nothing had ever called it from a `.cc` before).

   Standalone correctness test: `src/check/lj_crf_tile_kernel.t.cc`, one
   hand-built 32x32 self-tile from aladip's real atoms/positions/charges/
   LJ params (mask bit set for pairs in `[0.2nm, cutoff_long]`, excluding
   pathologically close would-be-excluded pairs whose huge r^-12
   repulsion amplifies float-vs-double rounding far beyond what's
   meaningful for testing kernel arithmetic), vacuum + rectangular, both
   passing within float precision (`tol = 1e-4`, appropriate for
   `FPL_TYPE = float` in the default mixed-precision build).

9. **Done (v1 scope: no perturbation/EDS -- multi-energy-group and virial
   landed later, see below).**
   `PLAN.md` §10 step 9: `interaction::CUDA_Nonbonded_Interaction`
   (`src/interaction/nonbonded/interaction/cuda_nonbonded_interaction.
   {h,cc}`), a real `Nonbonded_Interaction` subclass wired into
   `create_nonbonded.cc`'s `accelerator == gpu_cuda` branch (replacing the
   old, permanent mismatch of a real `CUDA_Pairlist_Algorithm` paired with
   the CPU `Default_Nonbonded_Interaction`, which silently ignored every
   real tile the pairlist built). `calculate_interactions()` drives
   `CUDA_Pairlist_Algorithm::prepare()`/`update()` then a new
   `compute_forces_energies()` passthrough (added to
   `CUDA_Pairlist_Algorithm`/`_Impl`) that runs `gpu::lj_crf_tile_kernel`
   over all four classified buckets (`solute_short/long`,
   `solvent_short/long`), accumulating into a persistent per-atom GPU
   force buffer and two energy accumulators (`m_iac`/`m_charge`/`m_force`/
   `m_e_lj`/`m_e_crf`, all built/sized once in
   `CUDA_Pairlist_Algorithm_Impl::init()` -- which is now actually called
   from `CUDA_Pairlist_Algorithm::init()`; it existed but nothing invoked
   it before), then adds (`+=`, matching `Forcefield::calculate_
   interactions()` zeroing `conf.current().force`/`energies` once up
   front, not per-`Interaction`) the result into `conf.current().force`
   and `energies.lj_energy[0][0]`/`crf_energy[0][0]`.

   (Originally also hard-errored on exactly one energy group and no
   virial -- both lifted later, see the entries after step 10 below.)
   Does not implement the CPU twin-range performance optimization of
   freezing long-range forces between pairlist rebuilds --
   `calculate_interactions()` fully reruns `prepare()`+`update()` (full
   candidate rebuild + reclassification) and recomputes *all* four tile
   buckets from current positions every call, regardless of
   `sim.param().pairlist.skip_step`. Numerically exact for any single
   evaluation (right after a rebuild, "recompute fresh" and "reuse what
   was frozen at this same rebuild" are the same numbers) so this doesn't
   compromise the correctness test below, but doesn't reproduce
   `skip_step`'s performance characteristic or (for a multi-step
   trajectory drifting across the `cutoff_short`/`cutoff_long` boundary
   between what would have been rebuild steps) its exact per-step
   energies either. Revisit before this runs production MD.

   Found and fixed two real bugs while wiring this, not just new work:
   - `check_forcefield.cc`'s finite-difference hessian check calls
     `Nonbonded_Interaction::calculate_interaction()` (singular)
     unconditionally on every `Nonbonded_Interaction`, including ones
     whose `init()` already hard-errored (aladip's own topology defines
     2 energy groups, tripping the gate above). The base implementation
     unconditionally indexes `m_nonbonded_set[0]`, guarded only by an
     `assert` that's compiled out under the release build's `-DNDEBUG` --
     since this class never populates `m_nonbonded_set`, that's a real
     segfault (reproduced via `aladip_cuda`, not hypothetical). Fixed by
     making `calculate_interaction` virtual (it wasn't) and overriding it
     to fail safely instead of indexing out of bounds.
   - Even with that fix, the same test still crashed: `aladip_cuda.t.cc`
     calls `Forcefield::init()` without checking its return value, so a
     hard-errored `init()` doesn't by itself stop `calculate_
     interactions()` from running afterward against `CUDA_Pairlist_
     Algorithm_Impl` state (`m_iac`/`m_charge`/`m_force`/etc.) that
     `init()` never allocated -- `cudaMemset` on a still-null,
     zero-capacity buffer. Fixed with an `m_initialized` guard: both
     overridden methods return a safe no-op instead of touching any GPU
     state if `init()` didn't complete. See `KNOWN_ISSUES.md`'s
     `aladip_cuda` entry for the full account.

   Correctness test: `src/check/cuda_nonbonded_interaction.t.cc` --
   unlike the step-6/step-8 tests (pairlist only; kernel only, hand-built
   tile), this drives the exact `CUDA_Pairlist_Algorithm` +
   `CUDA_Nonbonded_Interaction` pairing `create_nonbonded.cc` uses, and
   compares against `Standard_Pairlist_Algorithm`'s real pairlist summed
   through `Nonbonded_Term::lj_crf_interaction` directly (all four
   buckets), vacuum + rectangular, both passing within float precision.
   (Originally forced every atom into a single energy group to work around
   the v1 gate; now runs against aladip's real, unmodified 2-energy-group
   topology -- see the multi-energy-group entry after step 10 below.)

10. **Done (part 1 of 2): real twin-range cadence, matching CPU exactly.**
    `PLAN.md` §10 step 10 / §9.4's prerequisite: step 9 landed a
    deliberate simplification (`CUDA_Nonbonded_Interaction::calculate_
    interactions()` fully rebuilt + reclassified + recomputed *all* four
    tile buckets every single call, ignoring `skip_step` entirely).
    That's now fixed to match `nonbonded_set.cc`'s real semantics
    exactly: `pairlist_update = !(sim.steps() % skip_step)` (identical
    expression, same operand order) gates `CUDA_Pairlist_Algorithm::
    update()` (candidate rebuild + classification) and the long-range
    (`solute_long`/`solvent_long`) force/energy recompute; short-range
    (`solute_short`/`solvent_short`) still recomputes every call
    regardless. `CUDA_Pairlist_Algorithm_Impl::compute_forces_energies`
    gained a `bool recompute_long` parameter and three new members
    (`m_longrange_force`, `m_e_lj_long`, `m_e_crf_long`, sized once in
    `init()`) that hold the long-range contribution frozen between
    `pairlist_update` steps, added into every call's total unconditionally
    -- mirroring `nonbonded_set.cc`'s `m_storage.force +=
    m_longrange_storage.force` line for line. The zero-before-recompute
    for these three lives strictly inside the `if (recompute_long)`
    branch; zeroing them unconditionally would silently wipe the frozen
    value a non-rebuild step depends on.

    Existing single-evaluation tests (`pairlist_cuda_equivalence`,
    `lj_crf_tile_kernel`, `cuda_nonbonded_interaction`) are unaffected and
    still pass: at `sim.steps() == 0`, `pairlist_update` is always true
    regardless of `skip_step`, so a single `calculate_interactions()`
    call still does a full rebuild+recompute exactly as before. This
    piece doesn't need its own new test since it's only observable across
    multiple steps -- see step 10's second part below, which does drive a
    multi-step trajectory and would catch a regression here too.

    `skin` still doesn't do anything on the GPU path -- the candidate
    rebuild stays on the exact same cadence as classification (`M = N`
    in the design's own terms). Decoupling that is the second half of
    this step, not yet done as of this paragraph; see the next entry
    once it lands.

    **Done (part 2 of 2): skin-based candidate-rebuild decoupling.**
    `PLAN.md` §9.4's skin buffer drift test, closing the gap this
    section's part 1 flagged: the (expensive) candidate rebuild
    (`reorder()`+`build_candidates()`) is now decoupled from
    classification via a real displacement-tracked Verlet-buffer
    criterion, not just plumbed through as an inert field.

    `CUDA_Pairlist_Algorithm_Impl` gained `needs_candidate_rebuild()`
    (decides, given `sim.param().pairlist.skin`, whether a rebuild is
    also due this classification cycle) and `rebuild_candidates()`
    (wraps `reorder()`+`build_candidates()`, then snapshots the current
    GPU position mirror into a new `m_candidate_ref_pos` member so the
    next `needs_candidate_rebuild()` call has something to diff
    against). `CUDA_Pairlist_Algorithm::update()` calls
    `needs_candidate_rebuild()` first and only rebuilds if it says so,
    then always calls `classify_tiles()` regardless -- classification
    always re-evaluates exact current-position distances against
    whatever candidates currently exist, so a stale-but-still-safe
    candidate set produces identical classification output to a freshly
    rebuilt one.

    The Verlet criterion itself: `needs_candidate_rebuild()` returns
    true unconditionally if no rebuild has ever happened yet, or if
    `skin == 0.0` (an explicit early-return, not left as an emergent
    property of `2*0 >= 0`, so degenerate-skin behavior is obviously
    correct on inspection). Otherwise it launches a new kernel
    (`gpu::launch_max_displacement`, `gpu/cuda/interaction/nonbonded/
    kernels/displacement.{h,cu}`, same host/device split as
    `lj_crf_tiles.h` for the same reason -- `gpu::Periodicity` isn't
    host-compilable) that computes, per atom, the periodicity-wrapped
    (nearest-image, not raw-subtraction -- chargegroups get rewrapped
    into the box every step, so a raw subtraction would see a spurious
    ~one-box-dimension jump whenever an atom's chargegroup crosses a
    periodic boundary between the last rebuild and now) displacement
    since `m_candidate_ref_pos`, reduces to a single max via a small
    self-contained per-block-partial-max kernel + host-side linear scan
    (not `reduction.h`'s existing sum-only `calc_sum` -- that's a
    multi-launch-with-host-sync pattern, a needless cost paid at every
    classification-cadence check rather than just rebuild steps), and
    returns true iff `2 * max_displacement >= skin` -- the standard
    Verlet safety argument: two atoms could in the worst case have
    closed the gap between them by that much since the candidates were
    last built at `cutoff_long + skin`.

    Also added `candidate_rebuild_count()` (test-support only, on both
    `CUDA_Pairlist_Algorithm` and `_Impl`) so the drift test below can
    confirm skin is actually reducing rebuild frequency, not just
    numerically coinciding with always-correct behavior -- a
    `needs_candidate_rebuild()` bug that always returned true would
    still pass a pure force/energy comparison.

    Test: `src/check/cuda_skin_drift.t.cc` -- the first CUDA test in this
    tree to drive a real multi-step trajectory (every other test is a
    single evaluation). Two independently-constructed runs at the *same*
    `skip_step` (isolating skin's effect from part 1's already-covered
    skip_step cadence change), driven through an identical deterministic
    per-atom position-perturbation sequence: a `skin = 0.0` baseline
    (rebuilds every classification cycle) vs. `skin > 0` (decoupled).
    Since `classify_tiles()` always computes exact current-position
    distances regardless of candidate age, these must match EXACTLY at
    every step if skin's buffer is never exhausted -- confirmed for
    vacuum + rectangular, 30 steps, with the decoupled run skipping 4 of
    6 candidate rebuilds. Also confirms the `skin = 0.0` case rebuilds
    exactly every classification cycle (explicit degenerate check).

    A `skip_step = 1` run is deliberately *not* used as this test's
    reference (an earlier draft of this test did, and it failed --
    correctly: `skip_step` itself changes the physics via part 1's
    frozen long-range forces, so a `skip_step = 1` run and a
    `skip_step = 5` run are never expected to match at all, regardless
    of skin. That divergence was a real finding about this test's design,
    not a bug in the decoupling code -- worth remembering if this test
    is ever extended.).

Step 11 (re-porting `Leap_Frog_*<gpuBackend>`) is done -- see below.
Multi-energy-group and atomic-virial support are also done -- see below.
Perturbation/EDS is not started (a substantially larger piece of work: a
whole second pairlist/force code path, not a reduction shape change).

### Multi-energy-group support for `CUDA_Nonbonded_Interaction`

Lifts step 9's single-energy-group restriction, the blocker that used to
mask `aladip_cuda`'s virial gate entirely (aladip's real topology defines
2 energy groups, tripping the energy-group check before virial's own gate
was ever reached -- see `KNOWN_ISSUES.md`).

The tile kernel's energy reduction previously produced two flat scalars
(`e_lj_total[0]`, `e_crf_total[0]`) regardless of which atoms were
involved; the CPU reference (`nonbonded_innerloop.cc`) always buckets by
`[topo.atom_energy_group(i)][topo.atom_energy_group(j)]` into
`configuration::Energy::lj_energy`/`crf_energy`'s `[num_groups][num_groups]`
matrix. Closing that gap:

- `CUDA_Pairlist_Algorithm_Impl` gained `m_atom_energy_group` (per-atom
  energy-group index) and `m_num_energy_groups`, built once in `init()`
  next to the already-cached `m_iac`/`m_charge` (topology is static for a
  normal run) -- no changes needed to `gpu::Topology`/`TopologyView` at
  all, since `iac`/`charge` already bypass that mirror the same way.
- `m_e_lj`/`m_e_crf`/`m_e_lj_long`/`m_e_crf_long` became
  `num_groups^2`-sized (flattened `[gi * num_groups + gj]`), not scalars.
  `gpu::NbSimParams` gained `num_energy_groups`, passed by value into the
  kernel like every other constant.
- `lj_crf_tile_kernel` (`lj_crf_tiles.cu`) now buckets its reduction by
  energy-group pair instead of reducing straight to one scalar. Atoms
  within a tile are sorted by spatial cell, not original index -- energy
  groups are contiguous atom-index ranges, but that locality doesn't
  survive the sort, so there's no shortcut around a genuine per-thread
  group lookup (`atom_energy_group[a1]`/`[a2]`, indexed by original atom
  index exactly like `iac`/`charge` already are). Implemented via dynamic
  shared memory sized `2 * num_groups^2` `FPL_TYPE`s (LJ + CRF bucket
  matrices): zeroed at kernel entry, each active thread does one
  shared-memory `atomicAdd` into its `[eg_i][eg_j]` bucket, then one
  `atomicAdd` per bucket (not per pair) flushes to the global
  `e_lj_total`/`e_crf_total` arrays. Degenerates to exactly the old
  single-`atomicAdd`-per-tile cost when `num_groups == 1`; for larger
  group counts it trades the previous warp-shuffle reduction tree for
  shared-memory-atomic contention on however many buckets a tile
  actually touches -- a deliberate, small efficiency give-up for
  correctness at arbitrary group counts, consistent with this codebase's
  already-established per-pair-atomicAdd-for-forces precedent (see
  `CLAUDE.md`'s transfer-count-not-kernel-time philosophy). Not worth a
  fancier bucketed-shuffle scheme for the group counts GROMOS systems
  actually use (single digits).
- `CUDA_Pairlist_Algorithm_Impl::compute_forces_energies()` now writes
  directly into `conf.current().energies.lj_energy`/`crf_energy`'s real
  `[gi][gj]` matrix (same accumulation style as the CPU inner loop) instead
  of returning flat `double&` out-params threaded through three layers
  (`CUDA_Pairlist_Algorithm`'s wrapper, `CUDA_Nonbonded_Interaction::
  calculate_interactions()`) -- all three signatures simplified
  accordingly.
- The `init()` hard-error gate (`topo.energy_groups().size() != 1`) is
  gone; the virial gate (lifted separately, see below) and
  perturbation/EDS gate immediately after it were untouched at the time
  this landed.

Correctness test: `cuda_nonbonded_interaction.t.cc` no longer forces a
single energy group -- it now runs against aladip's real, unmodified
2-energy-group topology and compares every `[gi][gj]` entry of the
resulting matrix against the direct CPU (`Standard_Pairlist_Algorithm` +
`Nonbonded_Term::lj_crf_interaction`) reference, bucketed the same way.
Passes exactly, vacuum + rectangular. `lj_crf_tile_kernel.t.cc` (kernel-
only) updated for the new signature, kept at 1 group (unchanged intent).

### Atomic virial support for `CUDA_Nonbonded_Interaction`

Lifts the `sim.param().pcouple.virial != math::no_virial` hard-error gate
in `init()`. Once it's gone, `aladip_cuda`'s "molecular virial (finite
diff)" check moves from the virial error message to the perturbation
one -- confirming this was genuinely the last thing standing between that
check and perturbation's own, separate, not-yet-started gate.

The CPU reference (`nonbonded_innerloop.cc`) accumulates the *atomic*
virial unconditionally, every pair, regardless of whether a virial was
even requested: `storage.virial_tensor(b, a) += r(b) * force(a)` for
`a, b` in `0..2`, where `r` is the same periodicity-wrapped pair distance
already computed for the LJ/CRF math and `force(a) = f * r(a)`. Unlike
energy, this is *not* bucketed by energy group -- one flat 3x3 tensor
regardless of group count. Molecular virial (what the failing check
actually needs) is a *separate* correction applied afterward by
`Molecular_Virial_Interaction` (`create_forcefield.cc`, pushed into the
forcefield sequence whenever `pcouple.virial == math::molecular_virial`,
completely independent of which `Nonbonded_Interaction`/accelerator filled
in the atomic virial) -- so the only thing `CUDA_Nonbonded_Interaction`
itself needs to get right is the atomic sum; the molecular correction
applies to it transparently with zero CUDA-specific code.

Implementation, following exactly the shape of the energy-group work
above but simpler (no per-group bucketing -- always exactly 9 elements):
- `CUDA_Pairlist_Algorithm_Impl` gained `m_virial`/`m_virial_long`
  (9-element `double` accumulators, mirroring `m_e_lj`/`m_e_lj_long`'s
  short/long-range-frozen split exactly), sized once in `init()`.
- `lj_crf_tile_kernel` (`lj_crf_tiles.cu`) already computes `rvec` (the
  pair distance) and `fr` (the force vector) per active pair for the
  existing LJ/CRF math -- the virial add-on is 9 more shared-memory
  `atomicAdd`s per pair (`rvec.x*fr.x`, `rvec.x*fr.y`, ... unrolled rather
  than indexed, since `FPL3_TYPE` has no `operator[]`) into a fixed
  9-element region of the same dynamic shared-memory block already used
  for energy bucketing, then one global `atomicAdd` per element (9 per
  tile, not per pair) at the end -- same bucket-then-flush pattern, just
  a fixed size instead of `num_groups^2`.
- `compute_forces_energies()` adds `m_virial + m_virial_long` directly
  into `conf.current().virial_tensor(b, a)` (`+=`, matching
  `Forcefield::calculate_interactions()` zeroing it once up front, same
  as force/energy) at the end, right next to the energy-matrix
  accumulation loop.
- The `pcouple.virial != math::no_virial` gate in `init()` is gone.

Correctness tests: `lj_crf_tile_kernel.t.cc` (kernel-only) gained a
9-element virial check against a direct CPU reference built the exact
same way the existing force/energy reference already was in that file --
passes exactly, vacuum + rectangular. End-to-end: `aladip_cuda`'s
"molecular virial (finite diff)" check no longer fails on the virial gate
(confirmed via the actual error message, which is now the perturbation
one instead) -- full molecular-virial finite-difference validation
against a real, non-trivial topology, not just a hand-built kernel test.

### Step 11: re-porting `Leap_Frog_*<gpuBackend>` (PLAN.md §10 step 11)

Prerequisite: PLAN.md §3.2's GPU mirror cache (`sim.cuda().
configuration_view()`/`topology_view()`), landed in a prior commit on this
branch. The old `leap_frog_gpu.cu` predated that cache and called
now-removed APIs (`conf.m_gpu`, `conf.copy_pos_vel_to_gpu()`) -- it was
disabled (`is_supported_backend` excluded `gpuBackend`) and not wired into
any `CMakeLists.txt`. Replaced entirely, not patched.

- New `gpu::CudaManager::configuration_view(conf, sync_pos_vel, full_resync)`
  gained a third parameter: `full_resync = true` does a complete
  `copy_to_device()` (pos/vel/force/box/tensors) instead of the cheap
  pos+vel-only path. Needed because the integrator's velocity update reads
  `old().force`, which the existing cheap sync doesn't cover (force is
  otherwise only ever read fresh on the GPU side within the same nonbonded
  evaluation that produced it).
- New `gpu::Configuration::copy_pos_vel_from_device()` and
  `CudaManager::sync_configuration_from_device()`: the GPU-mirror-to-CPU
  direction, symmetric to the to-device path, for the one sync-back per
  step.
- New kernels, `gpu/cuda/algorithm/integration/leap_frog_kernels.{h,cu}`:
  plain per-atom elementwise updates (no boundary/periodicity dispatch
  needed here, unlike every pairlist kernel -- integration doesn't rewrap
  chargegroups).
- `algorithm/integration/leap_frog_gpu.cc` (renamed from `.cu` -- plain
  C++, no kernel syntax; `groalgorithm` has no CUDA language enabled, same
  constraint `remove_com_motion_gpu.cc` already works within) implements
  `Leap_Frog_Velocity<gpuBackend>`/`Leap_Frog_Position<gpuBackend>::apply()`,
  mirroring the CPU reference's exact math and state-exchange order
  (`leap_frog.cc`).

**The "no round trip between two algorithms" milestone, concretely:**
`Leap_Frog_Velocity<gpuBackend>::apply()` does the one real resync of the
step (`full_resync = true`, right after `conf.exchange_state()` so the
GPU mirror's `old()`/`current()` labeling matches the CPU's post-exchange
labeling exactly), writes the new velocity into the GPU mirror's
`current().vel`, and returns -- no sync back to the CPU.
`Leap_Frog_Position<gpuBackend>::apply()` then calls `configuration_view()`
with both sync flags `false`: it reads that same cached mirror (same
`conf.id()`, hits the 1-entry fast path), so it sees the velocity Velocity
just wrote without any re-upload, computes the new position on the GPU
mirror's `current().pos`, and only *then* calls
`sync_configuration_from_device()` once, since no other algorithm in the
sequence yet consumes GPU-resident state directly and `conf.current()`
must stay CPU-authoritative for everything downstream.

Polarisable (COS) charges are explicitly not supported on this path yet
(`Leap_Frog_Position<gpuBackend>::apply()` hard-errors via
`io::messages.add(..., io::message::error)` if
`sim.param().polarise.cos` is set, rather than silently computing the
CPU-only `posV` update wrong).

Test: `src/check/leap_frog_gpu.t.cc` -- two independent simulations built
from the same aladip topology/config, driven through 5 steps of
deterministic, externally-supplied per-atom forces (no real nonbonded
evaluation needed; this test is about the integrator, not the force
field), one run through the CPU backend and one through the GPU backend.
Positions and velocities must match within float tolerance at every step
(confirmed: exact match, 72 atoms, 5 steps).

### Follow-up: GPU-mirror freshness tracking (data-level cache coherence)

The "no round trip between two algorithms" design above shipped with the
sync decision at each `configuration_view()` call site hand-picked as a
`sync_pos_vel`/`full_resync` boolean -- each one encoding a private
assumption about what ran immediately before it. That produced a real,
silent-wrong-answer bug: `Leap_Frog_Position<gpuBackend>` assumed nothing
touches `current().vel` between it and `Leap_Frog_Velocity<gpuBackend>`, but
`create_md_sequence.cc` inserts CPU-only temperature coupling
(`Berendsen_Thermostat`/`NoseHoover_Thermostat`) between them whenever
`sim.param().multibath.couple` is set. That thermostat rescales
`conf.current().vel` on the CPU directly; the GPU-resident value Position
actually read never saw it. Wrong trajectory, no error message, whenever
temperature coupling is enabled -- the common case for real runs.

Fixed by replacing the boolean flags with data-level freshness tracking,
scoped to what's buildable without instrumenting every raw CPU write into
`configuration::Configuration` (that would mean wrapping `math::VArray`
itself, used by hundreds of CPU-only call sites and explicitly kept
GPU-unaware per this directory's `CLAUDE.md`):

- `src/gpu/mirror_fields.h` (new, zero CUDA dependency): `gpu::MirrorField`
  bits (`MIRROR_POS`/`VEL`/`FORCE`/`BOX`/`ALL`), granularity matching
  `ConfigurationStateView`, covering `current()`+`old()` together (every
  upload path already moves both in one call).
- `gpu::Configuration` (`configuration_struct.h`) gains two bitmask members:
  `gpu_fresh_fields` (this field's GPU copy is trustworthy, no resync
  needed for a read) and `gpu_dirty_fields` (subset of fresh -- the GPU
  holds a value the CPU hasn't seen yet, e.g. from `mark_gpu_dirty()`).
  Two different concepts: a field can be fresh (just uploaded, matches CPU)
  without being dirty (GPU-only, diverging from CPU).
- `CudaManager::configuration_view(conf, read_fields)` replaces the boolean
  signature: any requested field not already fresh gets resynced (the
  existing `copy_pos_vel_to_device()`/`copy_to_device()` granularities,
  chosen by tracked freshness instead of an assumption at the call site).
  `mark_gpu_dirty(conf, fields)` lets a kernel-writing caller (`Leap_Frog_
  Velocity<gpuBackend>`, for VEL) vouch for a field without a round trip.
- **Two new hooks around every `Algorithm::apply()` in `Algorithm_Sequence
  ::run()`** (also the special-cased `Remove_COM_Motion` call in `init()`):
  `flush_gpu_dirty(conf, touches)` BEFORE `apply()`, `invalidate_gpu_mirror
  (conf, touches)` AFTER. `touches` comes from a new virtual
  `Algorithm::gpu_mirror_touches()`, default `MIRROR_ALL` -- every ordinary
  (CPU-side) algorithm, including every not-yet-real "GPU-named" backend
  like `Remove_COM_Motion<gpuBackend>`, needs zero changes and is handled
  correctly by default. Two hooks, not one, because invalidating only
  *after* an algorithm ran isn't enough: if that algorithm reads/writes a
  field while it's still GPU-dirty, it operates on the stale pre-write CPU
  value *before* the after-hook ever fires. The before-hook publishes
  anything currently dirty in `touches` to CPU first, so the algorithm
  about to run sees the correct value. (This gap was caught by writing the
  regression test below -- the first attempt, with only an after-hook,
  reproduced a different, subtler wrong-answer bug: invalidating without
  first publishing made `Leap_Frog_Position<gpuBackend>` resync from a CPU
  copy that itself had never received Velocity's write, discarding it
  entirely instead of just missing the thermostat's contribution.)
- `Leap_Frog_Velocity<gpuBackend>` and `Leap_Frog_Position<gpuBackend>` both
  override `gpu_mirror_touches()` to return `0` -- neither ever writes
  `conf` directly (Velocity's `exchange_state()`/box copy aside, which is
  read fresh again next step regardless); everything goes through
  `sim.cuda()`'s tracked API. Returning `0` stops the after-hook from
  immediately erasing the VEL-fresh bit Velocity just set (which Position
  depends on reading without a resync), and stops the before-hook from
  wastefully flushing-then-re-uploading on Position's own call when nothing
  actually interleaved.

Regression test: `leap_frog_gpu.t.cc` gained
`run_thermostat_interleaved_case()` -- a real `algorithm::Algorithm_Sequence`
containing `[Leap_Frog_Velocity<gpuBackend>, a local `Vel_Scale_Algorithm`
(CPU-only, rescales `conf.current().vel`, default `gpu_mirror_touches()`,
exactly the shape of the real thermostat classes), Leap_Frog_Position
<gpuBackend>]`, run via `Algorithm_Sequence::run()` so the real hooks fire,
compared against the identical three-algorithm sequence on the CPU backend.
Before the fix: fails (velocity ratio off by a compounding `(1/scale_factor)`
factor per step -- the thermostat's contribution silently never reaches the
GPU-resident value Position reads). After: exact match, 3 steps, 72 atoms.
