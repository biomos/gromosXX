# Known issues

## `aladip_cuda` fails (by design, not a bug): perturbation isn't supported yet

- **Test affected:** `aladip_cuda` (`cuda-on` only), specifically the
  "molecular virial (finite diff)" check within it -- all other checks in
  that test pass. All other regression tests (`aladip`,
  `aladip_unperturbed`, `aladip_special`, `aladip_atomic`, `c16_cg`,
  `lambdas`, both presets) pass, as do the CUDA-specific
  `pairlist_cuda_equivalence`, `lj_crf_tile_kernel`,
  `cuda_nonbonded_interaction`, `cuda_skin_drift`, and `leap_frog_gpu`
  tests.
- **Why:** `CUDA_Nonbonded_Interaction::init()` hard-errors when
  `sim.param().perturbation.perturbation || sim.param().eds.eds` --
  this class never builds a `Perturbed_Nonbonded_Set` (or any
  `Nonbonded_Set` at all), and `CUDA_Pairlist_Algorithm::
  update_perturbed()` still errors and does nothing (no perturbed-
  pairlist support at all yet). `io::messages.add(..., io::message::
  error)`, so the run fails loudly with a clear message rather than
  quietly computing wrong or crashing. Confirmed by the actual error text
  observed running `aladip_cuda`: `"CUDA_Nonbonded_Interaction does not
  support perturbation or EDS yet"`.
- **Multi-energy-group and virial support both landed** (previously the
  blockers here, hit in that order before perturbation's gate was ever
  reached -- energy groups first, then virial, each unmasking the next):
  - `CUDA_Pairlist_Algorithm_Impl::compute_forces_energies()` writes
    directly into `configuration::Energy::lj_energy`/`crf_energy`'s real
    `[gi][gj]` matrix -- the tile kernel (`lj_crf_tiles.cu`) buckets its
    energy reduction by energy-group pair via a per-atom
    `atom_energy_group` array and dynamic shared-memory atomics
    (degenerating to the old single-atomicAdd-per-tile cost when there's
    only one group). Verified against aladip's real, unmodified
    2-energy-group topology (no synthetic override) by
    `cuda_nonbonded_interaction`, every `[gi][gj]` entry against a direct
    CPU (`Standard_Pairlist_Algorithm` + `Nonbonded_Term::
    lj_crf_interaction`) reference, vacuum + rectangular.
  - The same method also fills `conf.current().virial_tensor` with the
    atomic virial (`virial(b,a) += r(b)*force(a)`, exact CPU formula,
    not energy-group-bucketed -- a fixed 9-element shared-memory
    accumulator in the tile kernel, same bucket-then-flush pattern as the
    energy reduction). `Molecular_Virial_Interaction`
    (`create_forcefield.cc`) applies its center-of-mass correction to
    this generically, regardless of accelerator, with no CUDA-specific
    code needed for that part. Verified at the kernel level
    (`lj_crf_tile_kernel.t.cc`, all 9 tensor elements against a direct
    CPU reference, vacuum + rectangular) and end-to-end (this is what
    unmasked perturbation's gate above -- the virial-specific error
    message is gone, replaced by the perturbation one).
- **Found and fixed along the way, not just a design gap:**
  `check_forcefield.cc`'s finite-difference hessian check calls
  `Nonbonded_Interaction::calculate_interaction()` (singular) on every
  `Nonbonded_Interaction` unconditionally, including ones whose `init()`
  already hard-errored. The base implementation unconditionally indexes
  `m_nonbonded_set[0]` (guarded only by an `assert`, compiled out under
  the release build's `-DNDEBUG`) -- since `CUDA_Nonbonded_Interaction`
  never populates `m_nonbonded_set`, this was a real segfault, not a
  hypothetical one (reproduced while wiring §9, before the fixes below).
  Fixed two ways: (1) made `Nonbonded_Interaction::calculate_interaction`
  virtual (it wasn't) so `CUDA_Nonbonded_Interaction` can safely override
  it instead of crashing; (2) added an `m_initialized` guard to
  `CUDA_Nonbonded_Interaction` itself, since `aladip_cuda.t.cc` (like
  some other test harnesses) calls `Forcefield::init()` without checking
  its return value -- a hard-errored `init()` does not, by itself, stop
  `calculate_interactions()`/`calculate_interaction()` from being called
  afterward on `CUDA_Pairlist_Algorithm_Impl` state (`m_iac`/`m_charge`/
  `m_force`/etc.) that `init()` never allocated.
- **Not fixed here:** real perturbed-pairlist support requires
  `CUDA_Pairlist_Algorithm::update_perturbed()` (not started) and a
  perturbation-aware force path in `CUDA_Nonbonded_Interaction` (also not
  started) -- a substantially larger piece of work than energy groups or
  virial (a whole second pairlist/force code path, not just a reduction
  shape change). Until it lands, `aladip_cuda` is expected to fail. Don't
  try to make it pass by loosening the perturbation hard-error -- that
  would silently reintroduce exactly the "quiet fake result" failure mode
  `PAIRLIST_PLAN.md` §5(A) was written to avoid.
- Non-perturbed CUDA runs, with any number of energy groups and any
  virial setting, work end-to-end now (`TILE_PAIRLIST_DESIGN.md` §8/§9):
  real candidate-build + classification (steps 3-7, both chargegroup- and
  atomic-cutoff) feeding a real LJ + reaction-field force/energy/virial
  kernel with per-energy-group-pair bucketing (step 8), wired into a real
  `CUDA_Nonbonded_Interaction` that `create_nonbonded.cc` actually selects
  for `accelerator = cuda` (step 9, replacing the old always-CPU
  `Default_Nonbonded_Interaction` pairing).

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
