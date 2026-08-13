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
  block) work through the dummy pairlist today: they get a real, empty
  pairlist, a one-time warning, and zero nonbonded forces/energies -- also
  by design, also not yet physically meaningful, but at least not a crash
  and not silently wrong. There's no regression test for that path yet;
  worth adding one that specifically checks for the warning +
  zero-nonbonded-energy combination once `PLAN.md` §10 step 5 makes it
  worth asserting against (a passing/energy-bearing version of it will
  replace this note).
