# CUDA pairlist: findings and plan

Status: **findings written up, decision needed before implementation** (see
"Open decision" below). This document is scoped to one question:
*what should `create_nonbonded.cc` and a CUDA `Pairlist_Algorithm` look like
right now*, and *does the `Pairlist_Algorithm` interface need to change to
eventually host a real, tile-based GPU pairlist*. It complements `PLAN.md`
(the overall CUDA rewrite plan) rather than replacing any of its decisions;
where this document touches something `PLAN.md` already decided (D5, D6,
§6, §7), that decision is treated as fixed and just gets applied concretely.

## 1. How `create_nonbonded.cc` diverged from master, and why it crashes

`git diff master -- create_nonbonded.cc` shows exactly what changed. Two
things happened, both self-inflicted during CUDA prototyping, neither
intentional per the user:

1. The real pairlist-algorithm selection block was commented out:
   ```cpp
   // if (sim.param().pairlist.grid == 0) {
   //   pa = new Standard_Pairlist_Algorithm();
   // } else if (sim.param().pairlist.grid == 1) {
   //   pa = new Extended_Grid_Pairlist_Algorithm();
   // } else if (sim.param().pairlist.grid == 2) {
   //   pa = new Grid_Cell_Pairlist(topo, sim);
   // } else { ... error ... }
   ```
   replaced by an `io::messages.add("Ignoring pairlist setting, using CUDA
   pairlist algorithm", ...)` warning that doesn't construct anything. `pa`
   (declared as `Pairlist_Algorithm * pa = nullptr;`) is therefore **always
   null**, in every build, CPU or CUDA, regardless of `PAIRLIST`/`GRID` or
   `accelerator`. This is what segfaults all 7 regression tests (recorded in
   `KNOWN_ISSUES.md` from the previous cleanup pass — that entry is now
   superseded by this document's fix).
2. `pa->set_parameter(&ni->parameter());` was commented out, and
   `Pairlist_Algorithm::set_parameter()`/`m_param` (and the forward
   declaration `class Nonbonded_Parameter;`) were commented out in
   `pairlist_algorithm.h` to match. Checked: `m_param` is never actually
   *read* anywhere in any of the four CPU pairlist algorithms today — it's
   write-only. So this specific comment-out is inert (doesn't change
   behavior), but it's still a deviation from master with no purpose now,
   and the user asked for master's shape back.

Both were, per the user, a deliberate temporary hardcode to test something
during CUDA prototyping, not a design decision — so §4 below just restores
them.

## 2. What a `Pairlist_Algorithm` actually has to produce, and who reads it

`interaction::Pairlist` is `std::vector<std::vector<unsigned int>>`: a
per-atom row of partner-atom indices (`pairlist_solute[i]` = all `j` atom
`i` interacts with). `PairlistContainer` is six of these
(`solute_short/long/candidates`, `solvent_short/long/candidates`).

This shape is read *directly*, by index, in the innermost loop of
`Nonbonded_Outerloop::_lj_crf_outerloop` (and every other outerloop
variant — perturbed, QM/MM, electric field, latticesum):
```cpp
for (i = 0; i < end; ++i)
  for (j_it = pairlist_solute[i].begin(); j_it != pairlist_solute[i].end(); ++j_it)
    innerloop.lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity, ...);
```
This is a **scalar, host-memory, per-atom double loop**. `begin`/`end`/
`stride` (part of `Pairlist_Algorithm::update`'s signature) is how OMP
threads and MPI ranks split the atom range — the whole parallelization
model for the CPU nonbonded force calculation is built around being able to
slice this exact container by atom index.

Grep confirms `PairlistContainer`/`Pairlist_Algorithm` are used far beyond
just `Standard_Pairlist_Algorithm`: `mpi_nonbonded_master/slave`,
`perturbed_nonbonded_set`/`perturbed_nonbonded_outerloop`,
`qmmm_nonbonded_outerloop`/`qmmm_nonbonded_set`, `qm_zone.cc` all consume
this exact type. **This interface is the load-bearing wall of the entire
CPU nonbonded stack, not an implementation detail of one algorithm.**

## 3. What the GPU-native pairlist actually wants to produce

`gpu::TileContainerT<TileVecT>` (`gpu/cuda/memory/pairlist/tile.h`) mirrors
`PairlistContainer`'s six-way split structurally (same
`solute_short/long/candidates`, `solvent_short/long/candidates` names), but
each `TileVecT<Interaction_TileT<32,32>>` holds fixed-size 32×32 tiles: an
`index` (which 32-atom row-block × 32-atom col-block) plus a 32-bit-packed
exclusion `mask`. It lives in CUDA unified memory and is written by
`atomicAdd`-based `push_back` from device code.

This is not a cosmetic difference from `Pairlist`. A `Pairlist` row is
ragged, host-resident, and indexed per-atom; a `TileVecT` is fixed-block,
device-resident, and indexed per-32×32-block. Converting one to the other
means walking every tile, expanding each set bit back into an
`(i, j)` pair, and copying the whole thing host↔device — which is *exactly*
the transfer cost `PLAN.md` §1 identifies as the actual thing to avoid
("whole-step-on-GPU is the goal... redundant host-device transfers being
the big one"). A tile-native pairlist that has to materialize a
`PairlistContainer` to satisfy `Pairlist_Algorithm::update()`'s signature
has defeated its own purpose before a single force is computed.

**Third data point, for completeness:** the already-present (but unwired,
per `PLAN.md` §7) `nb_kernels.cu` expects a *third*, different shape again —
a flat `uint2* pairs` array, one thread per pair, not tiles at all. So as
of today there are three mutually-incompatible "pairlist" shapes sketched
across the tree (`Pairlist`, `TileContainer`, flat `uint2` pairs), and none
of them has a real producer wired to a real consumer yet. This document's
recommendation (§4) is partly about not adding a fourth.

## 4. Recommendation: don't redesign `Pairlist_Algorithm` — bypass it for the real GPU path

**Do not change `Pairlist_Algorithm`, `Pairlist`, or `PairlistContainer`.**
They're correct and heavily depended-upon for every CPU-executed nonbonded
path (plain, OMP, MPI, perturbed, QM/MM). Redesigning them to somehow also
describe tiles would either (a) genericize them into something too abstract
to stay simple/maintainable — directly against `PLAN.md` §1's "prefer
explicit, small, well-named classes" goal — or (b) bolt a tile-shaped escape
hatch onto a container whose entire existing contract is per-atom rows,
confusing every future reader of the four working CPU algorithms.

**Do introduce a second, separate, small interface for the tile-based GPU
pairlist, once it's built for real** — not a `Pairlist_Algorithm` subclass.
Something like:
```cpp
namespace interaction {
  class Cuda_Pairlist_Algorithm : public algorithm::Algorithm {
  public:
    virtual int build(topology::Topology & topo,
                       configuration::Configuration & conf,
                       simulation::Simulation & sim,
                       gpu::TileContainer & tiles) = 0;
  };
}
```
This is consumed directly by the real tile-based `CUDA_Nonbonded_Interaction`
(`PLAN.md` §7/§10 step 9), which is itself not a `Nonbonded_Interaction`
subclass either, for the same reason — it doesn't produce/consume
`PairlistContainer` or go through `Nonbonded_Set`/`Nonbonded_Outerloop` at
all. `create_nonbonded.cc`'s existing `#if defined(USE_CUDA)` branch already
special-cases which `Nonbonded_Interaction` subclass gets built; it's the
natural, already-present seam to also special-case which pairlist gets
built, with its own local variable — it does **not** need to keep assigning
into the generic `Pairlist_Algorithm * pa` for the real GPU path.

This means today's `interaction::CUDA_Pairlist_Algorithm<Backend>` (in
`cuda_pairlist_algorithm.h`, `: public Pairlist_Algorithm`) is not the
shape the *real* GPU pairlist should end up having. It's fine, and useful,
as the **dummy/placeholder** the user is asking for right now (§5) — it
lets `pa` be non-null and the whole existing CPU nonbonded machinery run
end-to-end with CUDA enabled, without touching `Nonbonded_Set`. But it
should be understood as temporary scaffolding, replaced outright (not
evolved) once `PLAN.md` §10 step 5 (port the real tile pairlist) lands, at
which point it's deleted and the `Cuda_Pairlist_Algorithm`/tile-native
`CUDA_Nonbonded_Interaction` pair takes over the `USE_CUDA` branch in
`create_nonbonded.cc` instead.

**So: no `Pairlist_Algorithm` redesign. Yes, a new (currently
nonexistent) interface for the real tile-based pairlist, added later,
alongside — not instead of — the existing one.**

## 5. Open decision: what should the *dummy* CUDA pairlist actually do?

Two honest options for what `CUDA_Pairlist_Algorithm<Backend>::update()`
does right now, before any real GPU pairlist exists:

- **(A) Explicitly empty.** `update()` clears the pairlist, emits one
  `io::messages` warning ("CUDA pairlist algorithm is a placeholder;
  producing zero nonbonded pairs until PLAN.md §10 step 5 lands"), and
  returns. Forces/energies from nonbonded interactions will be exactly
  zero whenever this path is selected. Cheap, obviously fake, impossible to
  mistake for a real result, good for testing that the *class
  hierarchy/dispatch/build* works (which is what "hardcoding" it was
  apparently for) without ever risking someone reading a plausible-looking
  but physically meaningless energy number.
- **(B) CPU-backed stand-in.** `update()` internally delegates to a real
  CPU algorithm (e.g. owns a `Standard_Pairlist_Algorithm` and forwards to
  it) so forces/energies come out physically correct, while the class
  still presents itself (name, selection path) as "the CUDA algorithm".
  Useful if the immediate goal is testing something *downstream* of the
  pairlist (integrators, energy conservation, output) with
  `accelerator=gpu_cuda` set, without caring that the pairlist itself isn't
  really on the GPU yet.

**Recommendation: (A).** It matches `PLAN.md` D5's spirit (no silent
CPU-pairlist-under-CUDA fallback becoming a permanent, load-bearing
default) even for a temporary placeholder, and it fails loudly (visibly
wrong energies + a warning) rather than quietly (plausible energies that
are secretly running on the CPU pairlist while everyone believes CUDA is
in use). If what's actually needed right now is (B)'s behavior for some
specific downstream test, that's better done by literally selecting a CPU
`Pairlist_Algorithm` for that test run (`accelerator=disabled` or an
explicit CPU-forcing test fixture) than by making the "CUDA" class quietly
CPU-backed.

**This is the one thing in this document that needs your decision before
I implement anything** — everything else below assumes (A).

## 6. Concrete plan (assuming (A) above)

1. **Restore `pairlist_algorithm.h`**: uncomment `class Nonbonded_Parameter;`,
   `set_parameter()`, `m_param` and its constructor initialization. No
   behavior change (confirmed unused elsewhere), just parity with master
   and un-blocking `create_nonbonded.cc`'s `pa->set_parameter(...)` call.
2. **Restore `create_nonbonded.cc`'s pairlist selection** to master's
   `grid == 0/1/2` block, uncomment `pa->set_parameter(&ni->parameter())`,
   remove the leftover `"Ignoring pairlist setting..."` warning. Wrap the
   restored block so it's only reached when *not* building the CUDA
   pairlist, i.e.:
   ```cpp
   if (sim.param().gpu.accelerator == simulation::gpu_cuda) {
   #ifdef USE_CUDA
     pa = new CUDA_Pairlist_Algorithm<util::gpuBackend>();
   #else
     io::messages.add("accelerator = cuda requested but this is a CPU-only build",
       "create_nonbonded", io::message::error);
     return 1;
   #endif
   } else {
     // the restored grid==0/1/2 block, verbatim from master
   }
   ```
   (The `#ifdef USE_CUDA` around the whole file's CUDA branch stays as it
   is now — this is just what goes *inside* it.)
3. **Implement `CUDA_Pairlist_Algorithm<Backend>::update()`/`prepare()`**
   as the explicit-empty dummy from §5(A): `pairlist.clear()` (already
   correctly sized/reserved by `Nonbonded_Set` before `update()` is
   called, same as every other algorithm), one `io::messages` warning
   (guarded so it fires once per run, not once per pairlist rebuild —
   e.g. a `bool` member set after the first call), return. Leave
   `update_perturbed()` as its current explicit "don't use this" warning
   (already reasonable), just switch its raw `std::cerr` to `io::messages`
   for consistency with the rest of the codebase.
4. **Leave `Default_Nonbonded_Interaction(pa)` as the CUDA branch's
   `Nonbonded_Interaction`, unchanged.** With an empty pairlist this just
   computes zero nonbonded forces/energies through the ordinary CPU
   outerloop — correct, if trivial, behavior for a container that's empty.
   This also means the CUDA branch stops needing its own special-cased `ni`
   construction *for now* (`PLAN.md` §7's TODO stays exactly where it is;
   nothing here changes what real tile-based `CUDA_Nonbonded_Interaction`
   needs to look like later).
5. **Rebuild and rerun `ctest` for both presets.** Confirm the SEGFAULTs
   in `KNOWN_ISSUES.md` are gone (they should be — `pa` is never null
   again) and that the CPU-pairlist-selected regression tests (`aladip`,
   `aladip_unperturbed`, etc.) now produce real results rather than
   crashing. `aladip_cuda` (the one CUDA-only test) will report zero
   nonbonded energy/force and a warning in its log — that's expected and
   correct for the dummy, not a new bug; note it in a
   `KNOWN_ISSUES.md` update (replacing the entry this document supersedes)
   so nobody mistakes it for a real regression later.
6. **Delete the now-resolved entry in `KNOWN_ISSUES.md`**, replacing it
   with a short note that `aladip_cuda` currently reports zero nonbonded
   energy by design (dummy pairlist), pointing at this document and
   `PLAN.md` §10 step 5 for when that stops being true.

None of the above touches `Pairlist_Algorithm`, `PairlistContainer`, or any
of the four working CPU pairlist algorithms. It also doesn't touch
`gpu::TileContainerT`/`Interaction_TileT` or attempt any real tile-pairlist
work — that's still `PLAN.md` §10 step 5, unchanged, and out of scope here.
