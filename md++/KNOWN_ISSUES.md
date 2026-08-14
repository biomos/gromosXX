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

## FIXED: CUDA context corruption after runtime `atomic_cutoff` toggle

- **Root cause found and fixed:** `CUDA_Pairlist_Algorithm_Impl::
  prepare_cog()` sized `m_cg_cog` to `topo.num_solute_chargegroups()`,
  but `prepare_cog_kernel` (`periodicity.cu`'s `prepare_chargegroup`)
  writes `cg_cog[cg_i]` for *every* chargegroup (`cg_i <
  topo.num_chargegroups()`, solute and solvent alike -- solute gets a
  true centre of geometry, solvent gets its first atom's position), and
  `classify_tiles()` reads `m_cg_cog` back for solvent candidates too,
  not just solute ones. Whenever any solvent chargegroups exist (true
  for essentially every real system, including aladip), this was an
  out-of-bounds `__global__` write on *every single* `prepare_cog()`
  call -- found via `compute-sanitizer --tool memcheck`, which reported
  it immediately on the very first call, with no `atomic_cutoff` toggle
  involved at all. The toggle wasn't causal; it's undefined behaviour
  that happens to sometimes land on memory that doesn't matter and
  sometimes doesn't, and toggling `atomic_cutoff` apparently perturbed
  allocation timing enough to make the "doesn't matter" case rarer in
  the original manual repro. Fixed by sizing `m_cg_cog` to
  `topo.num_chargegroups()` (matching what the kernel actually writes
  and what `classify_tiles()` actually reads), one line in
  `cuda_pairlist_algorithm_impl.cu`'s `prepare_cog()`.
- **How this was actually diagnosed:** the original manual repro
  couldn't be reproduced through `ctest` (`aladip_cuda.in` has
  `PERTURBATION` on, so it never reaches `CUDA_Pairlist_Algorithm::
  update()`/`prepare()` at all -- only the immediately-erroring
  `update_perturbed()`). Built a dedicated non-perturbed reproduction
  (`src/check/cuda_atomic_cutoff_toggle.t.cc`: loads
  `aladip_unperturbed.in`, forces `accelerator = gpu_cuda`, drives the
  real `Forcefield`/`CUDA_Nonbonded_Interaction` through the exact
  `atomic_cutoff` toggle sequence `check_atomic_cutoff` uses, repeated
  for 10 cycles) and ran it under `compute-sanitizer --tool memcheck`,
  which immediately pinpointed the exact kernel, exact buffer, and
  exact byte offsets of the invalid writes -- far more direct than
  guessing from the "everything fails afterward" symptom alone.
- **Verified:** `compute-sanitizer --tool memcheck` reports zero errors
  on `cuda_atomic_cutoff_toggle` (10 toggle cycles),
  `pairlist_cuda_equivalence`, `cuda_nonbonded_interaction`, and
  `cuda_skin_drift` after the fix (all previously had this same
  undetected OOB write happening silently on every `prepare_cog()`
  call -- none of them had visibly failed before, since the write
  happened not to corrupt anything that mattered on those particular
  runs, which is exactly the nature of this class of bug). Full `ctest`
  green on both presets except the unrelated, documented `aladip_cuda`
  perturbation gate. `cuda_atomic_cutoff_toggle.t.cc` is now a permanent
  regression test.
- **Worth remembering:** this bug predated and was unrelated to
  `atomic_cutoff` specifically -- it's worth periodically running the
  CUDA test suite under `compute-sanitizer --tool memcheck` even when
  nothing looks broken, since undefined-behaviour memory bugs like this
  one don't reliably announce themselves via `cudaGetLastError()` or a
  visibly wrong result.
