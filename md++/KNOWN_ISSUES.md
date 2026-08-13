# Known issues

## `aladip_cuda` fails (by design, not a bug): perturbation isn't supported by the dummy CUDA pairlist

- **Test affected:** `aladip_cuda` (`cuda-on` only). All other regression
  tests (`aladip`, `aladip_unperturbed`, `aladip_special`, `aladip_atomic`,
  `c16_cg`, `lambdas`, both presets) pass.
- **Why:** `aladip_cuda` runs the perturbed aladip system with
  `accelerator = cuda`, so it calls `Pairlist_Algorithm::update_perturbed`.
  `CUDA_Pairlist_Algorithm` is currently an explicit dummy placeholder (see
  `src/gpu/PAIRLIST_PLAN.md` §5/§6) that doesn't implement perturbed
  pairlists at all -- `update_perturbed` reports an `io::message::error`
  and does nothing, by design, so the run fails loudly rather than
  quietly computing wrong energies.
- **Not fixed here:** real perturbed-pairlist support requires the actual
  tile-based GPU pairlist (`PLAN.md` §10 step 5 onward), not the dummy.
  Until then, this test is expected to fail. Don't try to make it pass by
  loosening `update_perturbed`'s error -- that would silently reintroduce
  exactly the "quiet fake result" failure mode `PAIRLIST_PLAN.md` §5(A)
  was written to avoid.
- Non-perturbed CUDA runs (`accelerator = cuda` without a `PERTURBATION`
  block) exercise the real candidate-build + classification pipeline now
  (`TILE_PAIRLIST_DESIGN.md` §3 steps 1-5: chargegroup cog/cell build,
  Thrust sort-by-key into fixed 32-wide atom blocks, block bounding-sphere
  computation, block-pair candidate search, exclusion + short/long
  classification) before still ending in the same explicit dummy:
  `update()` clears the CPU-facing `PairlistContainer` and warns, so
  nonbonded forces/energies are zero, by design, not yet a crash and not
  yet physically meaningful either -- `m_tiles.solute_short`/`solute_long`/
  `solvent_short`/`solvent_long` are real, exclusion-checked, classified
  atom-pair tiles now, but nothing downstream (no force kernel) consumes
  them yet. Manually verified (not via `ctest` -- see next entry for why)
  that this runs to completion without error across many steps. There's
  no regression test for this path yet; worth adding one that checks for
  the warning + zero-nonbonded-energy combination, and separately for "did
  the pipeline actually run without a CUDA error," once `PLAN.md` §10's
  remaining work (a force kernel) gives it something real to assert
  against -- though the actual correctness gate is the pairlist-
  equivalence test (`TILE_PAIRLIST_DESIGN.md` §5), not this.

## Latent bug: CUDA context corruption after runtime `atomic_cutoff` toggle (not exercised by the current test suite)

- **What:** manually running `aladip_cuda` against a locally-modified,
  non-perturbed copy of `aladip_cuda.in` (to force-exercise the
  non-perturbed candidate-build path -- see previous entry) surfaced a
  real bug: after `check_forcefield.cc`'s "NonBonded: atomic cutoff"
  comparison check runs (which flips `sim.param().pairlist.atomic_cutoff`
  to `true` on an already-`init()`-ed `CUDA_Pairlist_Algorithm` and calls
  `update()`/`prepare()` again), **every subsequent CUDA call in the
  process fails** with `cudaErrorInvalidValue` ("invalid argument"),
  including calls with no relationship to `CUDA_Pairlist_Algorithm` at all
  (e.g. `gpu::Configuration::copy_to_device`). This looks like permanent
  corruption of the process's CUDA context/state, not a one-off failure at
  the toggle point itself.
- **What's already fixed, and isn't the root cause:** `CUDA_Pairlist_
  Algorithm::prepare()`/`update()` didn't originally re-check
  `atomic_cutoff` per call, only relying on `init()`'s one-time hard-error
  -- so a runtime toggle after `init()` would run `thrust::sort_by_key`
  against a chargegroup-mode cog/cell build that `prepare_cog()` had
  silently skipped that cycle. That's fixed (both methods now guard
  per-call, matching `prepare_cog()`'s existing behavior) and is a real,
  worthwhile fix on its own, but it only stops *this specific class* from
  running against mismatched state -- it does not explain or prevent the
  process-wide corruption, which persists and affects unrelated code
  (`copy_to_device`) too.
- **Why this isn't in `ctest` today:** the real, checked-in `aladip_cuda.in`
  has `PERTURBATION` enabled, and (separately confirmed) a perturbed run
  never reaches `CUDA_Pairlist_Algorithm::update()`/`prepare()` at all --
  only `update_perturbed()`, which errors out immediately (see the first
  entry above). So today's test suite cannot hit this even though the
  underlying corruption is real and reproducible.
- **Not fixed here:** the actual root cause (why does exercising
  atomic_cutoff, even now that it correctly no-ops in
  `CUDA_Pairlist_Algorithm`, corrupt the CUDA context for the rest of the
  process?) is not understood yet. Candidates worth checking first:
  whether `check_forcefield.cc` constructs a second, separate
  `Nonbonded_Interaction`/pairlist-owning object for the atomic-cutoff
  comparison (rather than reusing the existing one) and something about
  two such objects coexisting corrupts shared CUDA/Thrust/CUB state; or
  a CUDA version/driver quirk specific to this environment (CUDA 13.3).
  Needs a focused repro (a minimal standalone program constructing two
  `CUDA_Pairlist_Algorithm`-like objects) before guessing further.
