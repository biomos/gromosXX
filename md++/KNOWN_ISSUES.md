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
  classification), for both chargegroup-cutoff and atomic-cutoff
  (`sim.param().pairlist.atomic_cutoff`, `TILE_PAIRLIST_DESIGN.md` §4.2/
  step 7) -- verified against `Standard_Pairlist_Algorithm` bit-for-bit by
  the `pairlist_cuda_equivalence` test (vacuum + rectangular, both cutoff
  modes, 4 cases total). `update()` still ends in the explicit dummy for
  the CPU-facing `PairlistContainer` it's actually asked to fill (clears
  it and warns), since nothing downstream (no force kernel) consumes the
  real tiles (`m_tiles.solute_short`/etc) yet -- so nonbonded forces/
  energies through the normal `Nonbonded_Interaction` path are still
  zero, by design.

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
- **What's already fixed, and isn't the root cause:** as of
  `TILE_PAIRLIST_DESIGN.md` §4.2/step 7, `atomic_cutoff` is now a fully
  supported, real code path (not a no-op) -- `prepare_cog()`/`reorder()`/
  `build_candidates()`/`classify_tiles()` all re-read
  `sim.param().pairlist.atomic_cutoff` fresh on every call, so a runtime
  toggle between calls is no longer "running against mismatched state"
  the way it was when this bug was first observed (a no-op'd chargegroup
  build followed by a `thrust::sort_by_key` that assumed it had run).
  That class of problem is gone. It does not explain or prevent the
  process-wide corruption described below, which persists and affects
  unrelated code (`copy_to_device`) too -- this entry stays open.
- **Why this isn't in `ctest` today:** the real, checked-in `aladip_cuda.in`
  has `PERTURBATION` enabled, and (separately confirmed) a perturbed run
  never reaches `CUDA_Pairlist_Algorithm::update()`/`prepare()` at all --
  only `update_perturbed()`, which errors out immediately (see the first
  entry above), regardless of `atomic_cutoff`. `check_atomic_cutoff`'s
  toggle (`check_forcefield.cc`, exercised by `aladip_cuda.t.cc`) still
  only ever hits `update_perturbed()` on this system for the same reason,
  confirmed by re-running the full `aladip_cuda` ctest after implementing
  real `atomic_cutoff` support (step 7): same pre-existing failure mode
  (wrong energies from the always-empty perturbed pairlist), no crash, no
  new failure. So today's test suite still cannot hit this even though the
  underlying corruption is real and reproducible, and real `atomic_cutoff`
  support landing doesn't change that -- but a *future* non-perturbed CUDA
  regression test that exercises a runtime `atomic_cutoff` toggle (there
  isn't one yet) would be worth checking against this bug specifically,
  since `atomic_cutoff` is no longer a guaranteed no-op that could mask it.
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
