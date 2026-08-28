# SUMMARY.md — CUDA port: changes, correctness, performance

Covers `2e9a9f7` (start of the CUDA port effort) through `146556965`
(tip of `cuda_claude_energy_mirror`, current). 80 commits across the
topic branches listed below. Written for reporting purposes — precise
claims only, sourced from commit messages, `PERFORMANCE.md`,
`KNOWN_ISSUES.md`, and this session's own directly-measured benchmarks.
Numbers not yet independently retested under a consistent methodology
are marked **[pending retest]** — see "Performance retest plan" at the
end.

## Scope and methodology caveats

- **Not all 80 commits carry an independent performance measurement.**
  Many are feature-completeness or infrastructure steps (porting one
  more algorithm to GPU, refactors, bug fixes found by a new test) with
  no standalone perf claim attached at the time. This document reports
  correctness for each phase and performance only where it was actually
  measured, rather than inventing per-commit numbers.
- **Only one CPU baseline has been measured with OpenMP enabled**, and
  only starting from commit `f10578e6a` (this session). Everything
  benchmarked before that compared GPU against **single-threaded** CPU
  — not a fair comparison, and explicitly flagged as such below. A
  single, authoritative CPU+OMP reference run (master branch,
  `build_type=Release`, OMP on) is planned but **not yet executed** —
  see the retest plan.
- **The benchmark environment is noisy.** `nvcc` crashes intermittently
  (internal compiler assertion, unrelated to code correctness — retried
  until it succeeds). GPU run-to-run variance is real (~1-2%), and at
  least one test (`ubiquitin_gpu`) is provably sensitive to GPU
  contention from concurrently running processes (passed once under
  parallel `ctest -j4` by chance, then failed deterministically 3/3
  times standalone before a real fix; see the virial-bug entry below).
  All performance numbers in this document that come from this
  session were measured standalone, not under concurrent load; historical
  numbers from `PERFORMANCE.md` should be assumed to be standalone too
  but weren't re-verified for that condition.
- **Correctness rigor varies by commit.** The overwhelming majority of
  correctness claims below mean "the relevant `ctest` target(s) pass at
  a stated tolerance" (typically `1e-4`–`1e-6` relative, `5e-3` for
  bonded terms near a formula's near-cancellation point) plus
  `compute-sanitizer --tool memcheck` reporting zero errors. Exactly one
  item (the molecular-virial bug, below) was additionally verified
  against a real, independently-generated CPU reference number to ~4
  significant figures. This document does not claim uniform rigor
  across all 80 commits.
- **One test fails throughout this entire history and still does:**
  `aladip_cuda`'s "molecular virial (finite diff)" check, because
  `CUDA_Nonbonded_Interaction` does not support perturbation/EDS. This
  is a known, documented scope boundary (`KNOWN_ISSUES.md`), not a
  regression. Every other checkpoint's "ctest passes" claim below
  excludes this one, consistently.

## Known limitations (as of the current tip)

GPU-accelerated (`accelerator=gpu_cuda`) runs do **not** support:
perturbation / free-energy calculations, EDS, GAMD, MPI, `Crossdihedral`
forces (falls back to CPU, correctly, but no GPU port exists),
distance/dihedral/angle/J-value/RDC restraints and local elevation
(CPU-only, no GPU port), coarse-grained/shifted-CRF models, triclinic or
truncated-octahedral boundary for `Lattice_Shift_Tracker`/pairlist,
multiple GPUs. `Position_Restraint` (`POSITIONRES`) and standard
bond/angle/dihedral/improper/nonbonded/SHAKE/SETTLE/LINCS/M-SHAKE/
temperature/pressure coupling are supported and GPU-native.

## Phase 1 — CUDA build infrastructure and pairlist foundation

*Checkpoints: `cuda_claude` (`40e87a97a`) → `cuda_claude_cleanup`
(`274935fb4`) → `cuda_claude_tile_pairlist` (`e98b0c96e`)*

Established the CMake/Ninja CUDA build (both `USE_CUDA=ON/OFF`
presets), fixed early link/segfault issues, added the tile-based GPU
pairlist's core data structures (`TileVecT`, cell/block candidate
search, exclusion CSR, atomic-cutoff axis), and de-templated
`CUDA_Pairlist_Algorithm` (dead `Backend` axis removed — it's never
dispatched polymorphically, always constructed directly under
`accelerator=cuda`).

**Correctness:** `pairlist_cuda_equivalence` (new in this phase) is the
first real numeric validation of anything built so far — asserts the
GPU tile pairlist and the CPU reference pairlist produce identical
per-atom short/long pair sets. Getting it to pass found two real bugs:
(1) `TileVecT` passed by value into kernels caused a use-after-free
(temporary's destructor freed shared buffers while the kernel could
still be running) — fixed with a non-owning `View` type; (2) a capacity
estimate truncated instead of rounding up, silently under-allocating,
and the resulting overflow error was queued but never displayed because
the test never flushed `io::messages` — fixed both the truncation and
the test's message handling.

**Performance:** not measured — no force/energy kernel existed yet at
this phase (pairlist only).

## Phase 2 — Twin-range cadence and skin decoupling

*Checkpoint: `cuda_claude_skin_cadence` (`e4020b3f4`)*

Real GROMOS twin-range semantics: long-range force/pairlist frozen
between classification cycles (`skip_step`-gated), and candidate-list
rebuild further decoupled from classification via a Verlet-buffer
displacement check against `PAIRLIST`'s `SKIN` parameter (existed as a
no-op field before this).

**Correctness:** `cuda_skin_drift` (new) drives a real multi-step
trajectory (30 steps) and requires a skin-decoupled run (fewer
candidate rebuilds) to match a skin=0 (rebuild-every-cycle) run
*exactly*, since classification always recomputes exact distances
regardless of candidate age. Passed, 4/6 candidate rebuilds skipped
with no numeric difference.

**Performance:** not measured at this phase. (This same `SKIN`
mechanism was later found to be silently inactive in the actual
benchmark config and fixed for real — see `126e8beba` in Phase 8.)

## Phase 3 — First real GPU algorithms (vanilla, low-risk)

*Checkpoints: `cuda_claude_atomic_cutoff` (`f66a5155d`, `Remove_COM_Motion`)
→ `cuda_claude_temperature_thermostat` (`8ecd9351c`, `Temperature_Calculation`
+ `Berendsen_Thermostat`)*

First algorithms ported to genuine GPU kernels (not just a `gpuBackend`
label wrapping the CPU loop). Established the reduction-kernel pattern
(local grid-stride accumulate + one `atomicAdd` per output element)
reused by every later reduction-shaped port.

**Correctness:** each has a dedicated CPU-vs-GPU comparison test
(`remove_com_motion_gpu`, `temperature_gpu`), zero `compute-sanitizer`
errors, full `ctest` green (aladip_cuda exception, as always). One real
bug found and fixed: `topo.temperature_groups()`'s CSR convention
differs from `energy_groups()`'s (leading-0/exclusive-end vs.
inclusive-end/implicit-start-0) — caught by a failing test, not by
inspection.

**Performance:** not measured at this phase.

## Phase 4 — All four bonded force terms

*Checkpoints: `5dc8d6047` (QuarticBond) → `527515dc7` (Angle) →
`765bd90e6` (ImproperDihedral) → `cuda_claude_bonded_terms`/`fa04007d1`
(Dihedral)*

`CUDA_{QuarticBond,Angle,ImproperDihedral,Dihedral}_Interaction`: one
thread per term (not per atom — term lists are the sparse/static
input), direct kernel launch every step, no pairlist involved. All four
hard-error on perturbation/GAMD (out of scope, consistent gate
throughout the project).

**Correctness:** each has its own CPU-vs-GPU test
(`quartic_bond_gpu`/`angle_gpu`/`improper_dihedral_gpu`/`dihedral_gpu`),
zero sanitizer errors, full `ctest` green.

**Performance:** not measured at this phase (measured much later, after
several intervening sync-tax fixes — see Phases 8–10).

## Phase 5 — Constraint algorithms (SHAKE, SETTLE, LINCS)

*Checkpoints: `cuda_claude_shake`/`ae78638b4` (solvent SHAKE) →
`cuda_claude_solute_shake`/`fcc140b3a` (solute SHAKE) →
`cuda_claude_settle`/`a5ce7dbf9` (SETTLE) →
`cuda_claude_lincs`/`ac5c84e9e` (LINCS, solute+solvent)*

Three genuinely different algorithmic shapes: solvent SHAKE is
per-molecule-independent (one thread per molecule, same sequential
Gauss-Seidel as CPU, bit-comparable); solute SHAKE is Jacobi-style
(one thread per constraint per round, converges to the same manifold,
not bit-comparable to the CPU's in-place sweep); SETTLE is closed-form,
no iteration (bit-comparable); LINCS's CPU recursion is already
Jacobi-shaped, so its GPU port needs no reformulation (bit-comparable).

**Correctness:** each has its own test at an appropriate tolerance
(`shake_gpu`/`solute_shake_gpu` looser due to non-bit-comparable
convergence order; `settle_gpu`/`lincs_gpu` at `1e-6` relative). One
real bug fixed along the way: `CUDA_Shake` double-constraining solvent
when `NTCS != shake` (`7be2b7414`).

**Performance:** not measured at this phase.

## Phase 6 — Pressure/temperature coupling completion, cleanup, feature gating

*Checkpoints: `cuda_claude_pressure_barostat`/`ce3f98d91` →
`cuda_claude_nosehoover`/`6c72a795c` → `cuda_claude_review`/`c40e6ecd7`
→ `cuda_claude_feature_checker`/`4232e7a9c`*

`Pressure_Calculation`/`Berendsen_Barostat` needed no GPU kernel at
this point (trivial host math, already correctly interleaved via the
existing mirror-flush hooks) — verified, not just assumed, with a
dedicated end-to-end test. `NoseHoover_Thermostat` reused
`Berendsen_Thermostat`'s existing kernel (same per-atom scale formula).
Closes out every "vanilla" single-step MD algorithm on the original
roadmap.

The review pass (`c40e6ecd7`) found two real, fixable inefficiencies
across all the above (bulk `memcpy` instead of per-atom struct-rebuild
loops for double-precision constraints; sparse instead of
`O(num_atoms)` force-download for bonded terms/`CUDA_Shake`) and
de-duplicated copy-pasted reduction/scale code between the two
thermostats.

`4232e7a9c` fixed a **user-reported bug**: the GPU
feature-compatibility matrix (`io::check_features()`) still reflected
the old, pre-port CUDA implementation and hard-errored on features that
had since landed (bonded terms, LJ/CRF, COM removal, multi-energy-group,
LINCS, SETTLE, etc.) — unlocked the real, tested feature set; left
locked, with reasons, everything genuinely unsupported.

**Correctness:** all of the above ctest-green + sanitizer-clean per the
standard pattern. `feature_checker_gpu` (new) directly reproduces and
verifies the fix for the reported bug.

**Performance:** not measured at this phase.

## Phase 7 — Real-scale regression test; found 5 more bugs

*Checkpoint: `a6b5939f3`*

Added `extended_test/ubiquitin` (~22.7k atoms, real multi-atom
chargegroups) as a permanent CPU/GPU regression test — every prior GPU
test used `aladip` (72 atoms), too small to exercise several real
failure modes. `ubiquitin_cpu` (bit-for-bit against a checked-in
reference, every build) + `ubiquitin_gpu` (step-0 checked against the
same CPU reference; rest of trajectory checked only for physical
sanity, since GPU/CPU trajectories are expected to chaotically diverge
on a system this size).

Getting the GPU run working at all found **5 bugs invisible at
`aladip`'s scale**, including: `thrust`/`cub` crashing at large `n` in
this specific CUDA/driver combination (replaced with hand-written
kernels); `TileVecT::size()` returning an unclamped counter that could
exceed real buffer capacity, causing an out-of-bounds read; a candidate-
capacity density estimate that silently dropped real pairs at scale
(fixed with a detect-and-retry path). Full list in `KNOWN_ISSUES.md`
history.

**Performance:** not measured (correctness milestone). This benchmark
topology/config is what all later performance numbers in this document
use.

## Phase 8 — GPU-resident data flow (the "no sync unless a real consumer needs it" redesign)

*Key commits: `78a5145a9` (rf_excluded) → `e122d51c7` (bonded/nonbonded
force → mirror) → `206d45635`/`dd1c2c718` (event-based GPU-GPU
ordering) → `aa16ad852` (dropped-force bug, user-reported) →
`7a836ec88` (position restraints + mixed-CPU/GPU warning) →
`48a6decc3`/`25682dd72` (deferred sync: Temperature_Calculation,
Leap_Frog_Position) → `50e986cbd` (PERFORMANCE.md established) →
`3c7d13684` (Berendsen_Barostat/RemoveCOMMotion/Lattice_Shift_Tracker
GPU port) → `d8cb8d5a3` (FPH double-precision accumulation) →
`a7150c5f0`/`f135e49bc` (CUDA_Shake/CUDA_Settle onto the mirror) →
`7bd23657d` (virial_tensor on-device atomicAdd)*

This is the architectural core of the whole effort: every GPU-native
writer (`atomicAdd`s force/virial/constraint_force directly into a
shared GPU-resident mirror, publishing to CPU only when a real
downstream consumer needs it (`needs_fresh_cpu_force()`,
`gpu_mirror_touches()`, `is_gpu_native()`), instead of syncing back
every call. Codified as a hard rule in `src/gpu/CLAUDE.md`.

**Two real bugs found and fixed, both from a user-reported production
run, both silent-corruption-not-crash:**

1. **`aa16ad852`** — any CPU-only special-force interaction
   (`POSITIONRES` etc.) running after a GPU-native bonded/nonbonded
   term read a stale/zero CPU force array (the mirror was
   GPU-dirty but nothing told `Forcefield` to publish it first), and
   that contribution was later silently wiped out entirely when the
   mirror was eventually flushed (overwrite, not merge). Surfaced as
   slowly escalating LINCS rotation warnings over thousands of steps —
   never caught by the existing test suite (zero coverage of
   `interaction/special/*` under CUDA at the time). Fixed by flipping
   the default (`needs_fresh_cpu_force()` now defaults `true`; only
   known GPU-native writers opt out) and adding a regression test.
2. **`7a836ec88`** — same root cause's performance half:
   `POSITIONRES` was now *correct* but forced a full mirror flush every
   step for one CPU-only term. Closed by porting `Position_Restraint`
   to GPU-native, plus a new generalized diagnostic
   (`is_gpu_native()`/mixed-CPU/GPU warning at startup) so any *future*
   CPU-only algorithm under GPU acceleration is visible up front instead
   of discovered as an unexplained slowdown or, worse, silently.

**`d8cb8d5a3` is a real, deliberate correctness-over-speed tradeoff:**
switched force/constraint_force/virial/kinetic-energy/pressure
accumulation from single- to double-precision (`FPH`) at the
`atomicAdd` summation point specifically (per-term compute stays single
precision) — eliminates float-accumulation error compounding across
hundreds of contributions per atom. Cost: real (`NonBonded`'s own
compute sub-line ~54s→~81s/10000 steps at the time — double-`atomicAdd`
runs at roughly half single-precision throughput on this GPU class).
Kept — judged a correctness fix, not a style choice.

**Correctness:** all ctest-green + sanitizer-clean per standard.
`d8cb8d5a3`'s note: 13 force-touching tests re-verified after the
precision switch specifically.

**Performance [pending retest for pre-OMP-baseline checkpoints]:**
first `PERFORMANCE.md` snapshot (`50e986cbd`) recorded, comparing GPU
against **single-threaded** CPU only — not the fair baseline used from
here on; kept for historical reference only:

| | ns/day | Wall (10000 steps) |
|---|---|---|
| CPU, single-threaded | ~2.2–2.25 | ~769–783 s |
| GPU (this phase, various points) | ~4.5–4.6 (pre-fixes) | ~840–920s-class, not precisely pinned per-commit |

Individual algorithm wins recorded during this phase (own-timer, 10000
steps, single-thread-CPU-baseline era): `Leap_Frog_Position` 30.9s→0.08s,
`Lattice_Shift_Tracker` 0.79s(CPU)→0.10s, `Berendsen_Barostat`
0.28s→0.08s. These are real but not on the OMP-comparable baseline —
retained for the historical record only.

## Phase 9 — Fair CPU baseline; a severe self-inflicted regression found and fixed

*Key commits: `f10578e6a` (OMP enabled by default) → `b70f950dc`
(event-leak fix) → `2a1e702f8` (PERFORMANCE.md update)*

Enabled OpenMP by default in both `cuda-on`/`cuda-off` CMake presets —
the CPU baseline from here on is `OMP_NUM_THREADS=6` (12 logical cores
available), not single-threaded. Re-benchmarking under this change
surfaced a **severe regression**: GPU wall time 375–381s → **~1300s**
(~3.5× slower), unrelated to OpenMP itself.

**Root cause:** `field_producer_events` (per-mirror-field list of CUDA
events, used for GPU-GPU stream ordering, introduced in this same
history's Phase 8 — `206d45635`/`dd1c2c718`) grew unboundedly on
`MIRROR_POS`/`MIRROR_VEL`/`MIRROR_VIRIAL`. Earlier in the project,
`MIRROR_POS`/`MIRROR_VEL` freshness was invalidated implicitly by a
"missing field" resync path; narrowing `Forcefield`'s and every
constraint algorithm's `gpu_mirror_touches()` toward `0` (this same
history's own "no unnecessary push" work) removed that implicit
clearing point without providing a replacement one. Result: ~5000
never-cleared CUDA events accumulated per field over 1000 steps,
each future `cudaStreamWaitEvent()` pass paying O(steps) — an O(steps²)
total cost. `MIRROR_VIRIAL` had a smaller, independent instance of the
same leak (new `atomicAdd`-merge call sites from `7bd23657d` never
paired with a clearing call). **This bug was introduced within this
project's own history (Phase 8) and found/fixed within it (Phase 9)** —
not a pre-existing or externally-inherited defect.

**Fixed:** `clear_stale_producer_events()` calls added at the correct
"this field starts fresh again" points (`Lattice_Shift_Tracker`,
`Leap_Frog_Velocity` for POS/VEL; `zero_mirror_force()` for VIRIAL).

**Correctness:** ctest green + sanitizer clean, standard pattern.
Directly confirmed via temporary instrumentation showing the event
count bounded after the fix (was growing ~5–6/step, unboundedly,
before).

**Performance — real numbers, 10000-step ubiquitin benchmark, GPU:
RTX 50-series class card:**

| State | ns/day | Wall time |
|---|---|---|
| CPU, single-threaded | 2.2–2.25 | 769–783 s |
| **CPU, `OMP_NUM_THREADS=6`** | **8.96–8.97** | **194.6–194.7 s** |
| GPU, event-leak regression (bug present) | ~1.3–1.4 (est.) | ~1300 s |
| GPU, event-leak fixed | 4.56–4.60 | 375.7–380.8 s |

At this checkpoint, **GPU is ~1.93–1.96× slower than 6-thread CPU** —
the first honest, apples-to-apples GPU-vs-CPU comparison in the
project's history.

## Phase 10 — Pairlist `SKIN` misconfiguration fix

*Checkpoint: `126e8beba`*

Profiled `CUDA_Pairlist_Algorithm::update()` directly (host-side
`std::chrono`; `nsys`/`ncu` were not usable in this sandboxed
environment — corrupted output). Found the split between candidate
rebuild and classification was near-even (13.9s/14.6s over a 2000-step
diagnostic) but **every classification cycle also triggered a full
candidate rebuild** — the `SKIN` decoupling mechanism built in Phase 2
was silently inactive because the benchmark's own `.imd` never set the
optional `SKIN` column (defaults to `0.0`, which forces rebuild-always).
A benchmark-configuration gap, not a code gap — `SKIN`'s own
correctness test (`cuda_skin_drift`, Phase 2) already covered the
mechanism itself.

**Fixed:** set `SKIN 0.6` (nm) in the benchmark config. Zero code risk.

**Performance — same 10000-step benchmark:**

| | ns/day | Wall time | `NonBonded` | pairlist sub-bucket |
|---|---|---|---|---|
| Pre-`SKIN` fix | 4.56 | 380.8 s | 249.2 s | 150.1 s |
| **Post-`SKIN` fix** | **5.30** | **328.1 s** | **186.6 s** | **86.5 s** |

Gap to CPU+OMP narrowed from **1.96× slower → 1.68× slower**.

## Phase 11 — GPU-native molecular virial: correctness bug found and fixed

*Checkpoint: `c59f99d17`, this session*

User directly flagged a pressure/virial discrepancy between CPU and
GPU runs on their own independently-generated benchmark output. Root
cause, confirmed via direct comparison against the CPU reference: a
**real, two-part correctness bug**, present since virial accumulation
moved on-device (`7bd23657d`, Phase 8) but never caught by any existing
test (none exercised molecular-virial pressure coupling under GPU with
this specific interaction ordering):

1. `Molecular_Virial_Interaction` (CPU-only, corrects atomic virial →
   molecular virial) read `conf.current().virial_tensor` without ever
   being told to flush it from the GPU mirror first — operated on a
   stale/zero value.
2. Even after fixing the read, the correction was silently clobbered:
   `Pressure_Calculation`'s own later mirror flush (`copy_constraint_
   data_from_device()`) *overwrites* `virial_tensor` from the mirror
   rather than merging, and the mirror itself was never told about the
   CPU-side correction.

**Fixed at the source, not patched around:** ported the correction to
a genuine GPU-native `CUDA_Molecular_Virial_Interaction` (per-
pressure-group center-of-mass reduction + correction, both on-device,
`atomicAdd`-merged into the mirror like every other GPU-native writer)
and `Pressure_Calculation<gpuBackend>` (reads `virial_tensor` via the
standard mirror-read protocol instead of the overwrite-prone flush).

**Correctness — directly quantified, not just "ctest passes":**

| | Pressure (step 0, ubiquitin) | vs. CPU reference |
|---|---|---|
| CPU reference | 1.15675×10¹ | — |
| GPU, bug present | 1.42163×10³ | **~123× too large** |
| **GPU, fixed** | **1.15681×10¹** | **matches to ~4 significant figures** |

Full virial tensor was symmetric (~1.5×10⁵ magnitude) when buggy vs.
correctly asymmetric (~10⁴ magnitude, matching the expected effect of
the center-of-mass correction) when fixed. A secondary finding during
this investigation: `ubiquitin_gpu` passed once by chance under
parallel `ctest -j4` (GPU contention shifted a chaos-sensitive
trajectory) and failed deterministically 3/3 times standalone before
the real fix — the fixed version passes standalone, deterministically,
every time (verified 4/4 runs plus repeated full-suite runs).

**Performance:** no regression; not the goal of this fix. (Incidentally
measured net-neutral-to-slightly-positive on the 1000-step benchmark
used in Phases 12–14 below, since it eliminated one call site's
`cudaStreamSynchronize()`.)

## Phase 12 — Skip unnecessary force flush for empty `Crossdihedral`

*Checkpoint: `a6160431e`, this session*

`sim.param().force.crossdihedral` defaults to `1` with no way to
disable it via the `FORCE` block, so `Crossdihedral_Interaction` is
pushed into every `Forcefield` regardless of whether the topology
defines any crossdihedral terms (most, including this benchmark's,
don't). Its own compute is free when the term list is empty, but the
base class's `needs_fresh_cpu_force() == true` default forced a full
device-wide sync + force-array D2H copy before it, every step, for
nothing.

**Fixed:** cache "has terms" once at `init()` (topology is static for a
run), skip the flush when empty.

**Performance — 1000-step ubiquitin benchmark, standalone, 2 runs:**

| | Wall time | ns/day |
|---|---|---|
| Before | 33.1–33.3 s | 5.51–5.53 |
| **After** | **29.6–29.7 s** | **6.20–6.22** |

**~11% faster.** Also reduced downstream `NonBonded` contention
(20.7s→18.25s over 1000 steps) — one fewer forced device-wide sync
point in the step.

## Phase 13 — GPU-native energy accumulation, bonded/special terms

*Checkpoint: `d8f2baed7`, this session*

Profiling (CUDA-event kernel timing + wall-clock brackets around every
stage of `CUDA_Angle_Interaction::calculate_interactions()`) found
~480ms/1000 steps hiding in what looked like a trivial two-double host
merge loop — not kernel execution (~205ms measured via events) or the
explicit sync (~190ms measured via wall clock), but a genuine CUDA
Unified Memory page-fault migration: `gpu::cuvector` (every bonded/
special term's private per-energy-group scratch buffer) is
`cudaMallocManaged`, and a host read immediately after a GPU kernel
write faults and migrates the page every step, for as little as 2
doubles. Cost concentrated on whichever GPU-native term ran first each
step (CUDA's unified-memory fault handling does batched/global
bookkeeping on first post-kernel host touch), not evenly distributed.

**Fixed:** new shared, plain-device (never host-touched except by one
explicit `cudaMemcpy`) energy buffers on the GPU mirror
(`gpu::MIRROR_ENERGY`), mirroring `virial_tensor`'s existing
`atomicAdd`-merge design. `Angle`/`Dihedral`/`ImproperDihedral`/
`QuarticBond`/`Position_Restraint` now merge on-device and do zero host
synchronization in `calculate_interactions()`. `Energy_Calculation`
(`gpu_mirror_touches()` now `MIRROR_ENERGY`, was `0u`) does the single
real publish per step.

**Performance — 1000-step ubiquitin benchmark, standalone, 2 runs:**

| Term | Before | After |
|---|---|---|
| `Angle` | 690 ms | **17 ms** |
| `Dihedral` | ~1000 ms | **14 ms** |
| `ImproperDihedral` | ~250 ms | **14 ms** |
| **Wall time total** | 29.6–29.7 s | **27.8–27.9 s** |
| Performance | 6.20–6.22 ns/day | **6.57 ns/day** |

**~6% faster** on top of Phase 12 (~16% cumulative from the Phase-12
baseline).

## Phase 14 — GPU-native energy accumulation, nonbonded LJ/CRF

*Checkpoint: `146556965`, current tip*

Same pattern applied to `CUDA_Nonbonded_Interaction`'s
`lj_energy`/`crf_energy` (per-energy-group-*pair* matrix, both
short-range and frozen long-range contributions), which had the
identical managed-memory host-readback pattern at much larger scale —
this is the dominant nonbonded cost, not a peripheral one.

**Performance — 1000-step ubiquitin benchmark, standalone, 2 runs:**

| | Before | After |
|---|---|---|
| `NonBonded` "compute forces energies" sub-timer | ~12.3 s | **0.97 s** |
| **Wall time total** | 27.8–27.9 s | **19.3–19.4 s** |
| Performance | 6.57 ns/day | **9.67–9.71 ns/day** |

**~31% faster** on top of Phase 13. `Leap_Frog_Velocity`'s own timer
absorbed some of the freed GPU contention (2.5s→5.1s — same
"whichever algorithm syncs next inherits the wait" effect observed
elsewhere in this project), consistent with, not contradicting, the net
wall-time win.

**Cumulative, Phases 12–14 (1000-step benchmark, standalone):
33.1s → 19.3s, ~1.7× faster.**

**Correctness, Phases 13–14:** full `ctest` green (aladip_cuda
exception only) both phases; `compute-sanitizer` zero errors on every
touched test (`angle_gpu`, `dihedral_gpu`, `improper_dihedral_gpu`,
`quartic_bond_gpu`, `position_restraint_gpu`, `cuda_nonbonded_
interaction`, `cuda_skin_drift`, `rf_excluded_gpu`, `pairlist_cuda_
equivalence`, `lj_crf_tile_kernel`, `ubiquitin_gpu`). Five standalone
test files needed updating for the new publish point (`gpu::MIRROR_
ENERGY` flush + read from `conf.old()` instead of `.current()`,
matching `Pressure_Calculation`'s existing convention) — a real,
expected consequence of moving the publish point, not a correctness
regression; each was verified to fail before the corresponding source
fix and pass after.

## Negative result: `TileVecT::m_size`/`m_overflow` → `cudaMalloc` (not committed)

Explored whether the tile pairlist's `m_size`/`m_overflow` counters
(currently `cudaMallocManaged`, host-read every classification cycle)
had the same page-fault-migration cost as the energy buffers above.
Implemented the `cudaMalloc` conversion, verified correct
(`ctest`+sanitizer clean), benchmarked with an isolated 3-way
comparison (original / memset-removed-only / memset-removed+
`cudaMalloc`) — **no measurable performance difference** (all three
within ~33.4–34.1s noise band on the 1000-step benchmark). Consistent
with the page-fault cost scaling with touch *frequency* × data volume:
these counters are touched ~200 times/1000 steps (once per
classification cycle, not every step) versus energy buffers' every-step
touch. Kept in the working tree for a possible correctness-motivated
follow-up, **not committed, not counted in any performance claim above.**

## Performance retest plan

**Goal:** replace the mixed-provenance numbers above (different step
counts, different CPU-baseline eras, some single-run) with a
consistent, defensible set for the presentation.

**CPU reference (once):** build `origin/cuda` (or the commit closest to
where OMP-fair comparison starts being meaningful,
i.e. no earlier than the topic-branch tips through Phase 6) at
`CMAKE_BUILD_TYPE=Release`, `OMP=ON`, `USE_CUDA=OFF`. Run the
10000-step ubiquitin benchmark once (matches `PERFORMANCE.md`'s
existing methodology), `OMP_NUM_THREADS=6`. This is the single CPU
number every GPU checkpoint below is compared against.

**GPU checkpoints to retest** (branch tip or commit, each rebuilt clean,
3 standalone repeats, mean reported):

| # | Checkpoint | Commit | Why |
|---|---|---|---|
| 1 | `cuda_claude_feature_checker` | `4232e7a9c` | End of "all vanilla algorithms ported" phase, pre-mirror-redesign |
| 2 | pre-event-leak-fix | `dd1c2c718` or `aa16ad852` | Documents the regression magnitude directly (expect ~1300s) |
| 3 | post-event-leak-fix, pre-OMP-baseline-note | `b70f950dc` | First honest GPU number after the regression fix |
| 4 | post-`SKIN` fix | `126e8beba` | Last-known-good before this session's correctness/perf work |
| 5 | post-virial-fix | `c59f99d17` | Correctness fix; confirm no perf regression at 10000 steps |
| 6 | post-`Crossdihedral` fix | `a6160431e` | Phase 12 |
| 7 | post-bonded-energy-mirror | `d8f2baed7` | Phase 13 |
| 8 | current tip | `146556965` | Phase 14, full cumulative number |

Checkpoints 1–3 do not need the OMP enabled preset (`f10578e6a`)
There is only one CPU ground truth, the master branch. Compile the master branch with build_type=Release and OMP enabled. Run benchmark on this with OMP_NUM_THREADS=6. This will be the baseline against which we compare everything. The GPU enabled runs do not need the OMP_NUM_THREADS.

**Methodology for every checkpoint:** clean rebuild (`rm -rf build/cuda-on
&& cmake --preset cuda-on && cmake --build`), retry on `nvcc` internal-
compiler crashes (environmental, unrelated to code), run the 2000-step
ubiquitin benchmark 3× standalone (not concurrent with any other GPU
process), report mean wall time / ns/day, note stdev. Re-run `ctest` at
each checkpoint and record pass/fail (expect `aladip_cuda` as the only
exception throughout).

Wrap up all results in a single jupyter notebook GPU_summary.ipynb. Draw relevant plots comparing performance. Be clear, consise and comprehensible.
