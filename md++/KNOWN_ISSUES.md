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

## FIXED: real-scale (extended_test/ubiquitin) GPU run was silently wrong -- a cascade of five bugs, none caught by any small-topology test

Adding `extended_test/`'s real ubiquitin-in-water simulation (~22.7k
atoms, 775 solute + 7044 SPC solvent molecules, real multi-atom
chargegroups) as a GPU-vs-CPU regression test surfaced five distinct,
previously-undetected bugs -- every existing GPU test topology (aladip:
12 solute atoms, 60 solvent atoms) was too small to trigger any of them.
All fixed; `ubiquitin_cpu`/`ubiquitin_gpu` (see below) are now permanent
regression tests guarding against a repeat.

1. **thrust/cub crash at real scale.** `CUDA_Pairlist_Algorithm_Impl::
   reorder()`'s `thrust::sequence`/`thrust::sort_by_key` calls crashed
   with `cudaErrorInvalidDevice`/`cudaErrorInvalidValue` originating
   inside cub's internal `PtxVersion()`/`CurrentDevice()` device query
   (`cub/util_device.cuh`), reproducible on this environment's GPU
   (very recent architecture) + CUDA 13.3 combination, but *only* for
   large `n` (never triggered by any small topology). A minimal
   standalone repro (same flags, same problem size, even with
   `-rdc=true`) never reproduced it -- root cause not fully isolated,
   suspected to be an RDC multi-translation-unit device-link interaction.
   Fixed by replacing both calls with hand-written kernels: a trivial
   `sequence_kernel` and a stable bitonic sort (`bitonic_sort_by_key`,
   `block_pairlist.h`/`.cu`) -- eliminates the thrust/cub dependency
   from this call path entirely, sidestepping the issue regardless of
   its exact cause.
2. **`TileVecT::size()` returned an unsafe, unclamped count.**
   `push_back()`'s `atomicAdd` on the size counter happens *before* the
   capacity check, so on overflow the counter ends up larger than the
   actual buffer -- `size()` returned that raw counter, and a caller
   trusting it as "safe to index up to this" (exactly what
   `classify_tiles()`'s kernel-launch grid dimension does) read past
   the real allocation. Fixed by clamping both `size()` accessors
   (`TileVecT` and `TileVecT::View`) to `min(*m_size, m_capacity)`
   (`tile.h`). `was_overflowed()` still reports the raw condition
   unaffected.
3. **Candidate-tile capacity silently underestimated, dropping real
   pairs.** `estimate_candidate_capacity()`'s density heuristic (a
   sizing *estimate*, not a correctness mechanism -- `was_overflowed()`
   is the real safeguard) underestimated badly enough at ubiquitin's
   scale to trigger real, repeated overflow -- candidates silently
   dropped, the pairlist genuinely wrong. The overflow *was* being
   detected and reported via `io::messages`, but nothing had displayed
   the message queue yet when the run diverged, so it looked like
   silent corruption. Fixed with a retry: on detected overflow,
   re-reserve at the true worst-case capacity (every possible block
   pair -- always mathematically sufficient, so at most one retry is
   ever needed) and re-run the candidate search
   (`_build_candidates()`, `cuda_pairlist_algorithm_impl.cu`).
4. **`Temperature_Calculation<gpuBackend>`'s reduction kernel used
   O(num_groups) *dynamic shared memory*.** Fine for a handful of
   user-chosen energy groups, but "temperature groups" are ~one per
   solvent molecule for a real system -- ubiquitin's ~7045 groups needs
   275KB, blowing past the ~48KB default dynamic shared memory limit
   and failing the kernel launch itself (`cudaErrorInvalidValue`)
   before a single thread ran. Fixed by switching
   `group_velocity_reduce_kernel` to direct global atomics (no
   per-block shared bucket) -- no ceiling, and per-group contention
   stays low regardless of `num_groups` since each group is only a few
   atoms (`temperature_kernels.cu`).
5. **1,4-pair ("LJ exception") interactions were never implemented on
   GPU at all.** `topo.all_exclusion(i)` (used for tile-classification
   masking) is `exclusion(i)` UNION `one_four_pair(i)` -- so 1,4 pairs
   were correctly *excluded* from the regular tile-driven LJ/CRF sum on
   both CPU and GPU, but the CPU path adds them back in with their own
   scaled `cs6`/`cs12` parameters and Coulomb-scaled CRF
   (`Nonbonded_Outerloop::one_four_outerloop`, called unconditionally)
   -- nothing on the GPU side ever did this. Aladip's tiny test
   molecule has few/no 1,4 pairs contributing meaningfully; ubiquitin's
   real backbone has many, producing a large (~2.6x), previously
   unexplained solute-solute CRF mismatch that survived ruling out the
   pairlist (a dedicated CPU-vs-GPU pair-set comparison found zero
   mismatches), floating-point precision (bit-identical wrong answer
   at `FP_PRECISION=1` and a full-double `=3` rebuild), long-range
   caching, and the shared RF constants. Fixed by adding
   `gpu::launch_one_four` (`one_four_kernels.h`/`.cu`, same CSR-based
   pattern as `launch_rf_excluded`) and extending `gpu::LJParams`/
   `LJParamView` to carry `cs6`/`cs12` alongside `c6`/`c12`
   (`cuda_lj_params.h`/`.cu`).
   - **A sixth, genuinely separate bug found while validating #5:**
     `interaction::Nonbonded_Parameter::m_coulomb_scaling` had *no*
     default member initializer and was silently dropped by the copy
     constructor -- any `Nonbonded_Parameter` constructed without an
     explicit `set_coulomb_scaling()` call (production always calls it
     via `create_nonbonded.cc`; several GPU check tests that construct
     `CUDA_Nonbonded_Interaction` directly do not) read uninitialized
     memory. This produced two extremely confusing, build/memory-
     layout-dependent test failures after adding the 1,4-pair kernel
     (`cuda_nonbonded_interaction` failing for `rectangular` but not
     `vacuum`; `cuda_skin_drift` failing for `vacuum` but not
     `rectangular` -- exactly backwards from each other, and not
     reproducible via a clean rebuild, which ruled out stale build
     artifacts). Fixed with a `= 1.0` default (the standard, non-AMBER
     case) and a corrected copy constructor (`nonbonded_parameter.h`).

**Verified:** full `ctest` green on both presets except the unrelated,
documented `aladip_cuda` perturbation gate above (plus `rf_excluded_gpu`'s
pre-existing intermittent flake, next section -- confirmed via
`git stash` on an unmodified checkout that it predates all of this).
`compute-sanitizer --tool memcheck`: zero errors on the affected tests.
The real 100-step ubiquitin GPU run that originally motivated this
investigation now completes successfully with matching (within
tolerance) step-0 energies against the CPU reference.

**A seventh bug found immediately after, the hard way:** the initial
version of `ubiquitin_cpu.t.cc`/`ubiquitin_gpu.t.cc` only ever compared
energies read directly from `configuration::Configuration` in memory --
never the actual `.tre`/`.trc` bytes a real run writes to disk. Prompted
to check those too, a manual diff of real `program/md` output initially
looked fine, but wiring an equivalent check directly into the test
(`extended_test/tre_parser.h`/`.cc`, driving a real
`io::Out_Configuration` exactly like `program/md`'s own main loop) found
that `ubiquitin_runner.cc`'s in-memory collection was reading
`conf.current().energies`, which is all-zero by the time a step's
`Algorithm_Sequence::run()` returns -- `io::Out_Configuration::print()`
itself reads `conf.old().energies` (confirmed directly in
`out_configuration.cc`), since leap-frog's own old/current rotation has
already moved the just-computed values there. A now-abandoned earlier
workaround (run once and discard before the real loop) happened to
paper over this by accident, for reasons never fully explained, and
would have kept doing so silently forever had the file contents not
been cross-checked against memory. Fixed by reading `conf.old()`
instead (`ubiquitin_runner.cc`); `ubiquitin_cpu.t.cc`/
`ubiquitin_gpu.t.cc` now also parse the real `.tre` file and require
every written step to match the in-memory values exactly, plus scan the
`.trc` file's positions for NaN/Inf/gross corruption -- so a future
regression in `Out_Configuration`'s own formatting/precision, or in
which configuration object anything reads from, can't hide from this
test the way it hid from the original, memory-only version.

## Known, unresolved: `rf_excluded_gpu` intermittent (~1-in-10) flake

- **Symptom:** a single-atom force/energy mismatch, small in magnitude,
  appearing nondeterministically (~1 run in 10-15) when
  `rf_excluded_gpu` is run repeatedly back-to-back; never reproduces
  twice with the same input in a row otherwise, and `compute-sanitizer
  --tool memcheck` reports zero errors on this test.
- **Confirmed pre-existing, not introduced by any change in this
  session:** reproduced identically via `git stash` on a checkout from
  *before* the thrust/cub-crash and bitonic-sort work (see the section
  above) -- same failure rate, same character, on the original
  `thrust::sort_by_key`-based code.
- **Not yet root-caused.** A plausible but unconfirmed hypothesis: some
  genuine GPU-side floating-point non-associativity (atomic accumulation
  order) sensitive to scheduling jitter between runs, though this
  wouldn't fully explain why it's specifically this test and not others
  of similar shape. Low priority given the magnitude is small and no
  other test has ever shown a similar pattern, but worth investigating
  properly before it's mistaken for noise on some future, larger change.
