# Known issues

## All `src/check` regression tests SEGFAULT (both `cuda-on` and `cuda-off`)

- **Tests affected:** `aladip`, `aladip_unperturbed`, `aladip_special`,
  `aladip_atomic`, `c16_cg`, `lambdas` (both presets), plus `aladip_cuda`
  (`cuda-on` only).
- **Root cause:** `src/interaction/nonbonded/create_nonbonded.cc`'s
  `Pairlist_Algorithm * pa = nullptr;` is never assigned a real algorithm.
  The block that would set it (`if (sim.param().pairlist.grid == 0) { pa =
  new Standard_Pairlist_Algorithm(); } ...`) is entirely commented out
  under a `//// TEMPORARILY OFF ////` marker, replaced by a
  `io::messages.add("Ignoring pairlist setting, using CUDA pairlist
  algorithm", ...)` warning that doesn't actually construct anything.
  `pa` therefore stays `nullptr` and is later dereferenced, crashing
  every run that builds a `Nonbonded_Interaction`.
- **Introduced by:** commit `228cd1dae` ("First attempt to implement
  separate cuda nonbonded interaction"), predates the `cuda_claude_cleanup`
  branch. Not a regression from that branch's compile-fix work.
- **Not fixed here:** restoring pairlist-algorithm selection is
  intertwined with `src/gpu/PLAN.md` §6's pairlist redesign (which CPU
  algorithm gets chosen when, and how that interacts with the future
  CUDA-only pairlist gate). Out of scope for a compile-only cleanup pass;
  see `PLAN.md` §10 roadmap step 5 onward.
