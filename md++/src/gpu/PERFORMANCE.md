# PERFORMANCE.md — CUDA port performance state

Snapshot as of commit `25682dd72` (branch `cuda_claude_rf_excluded`).
Written after a benchmarking/optimization pass covering the async
stream-ordered redesign, the `POSITIONRES`-dropped-force bug fix, the
mixed-CPU/GPU warning mechanism, the `Temperature_Calculation`
deferred-sync fix, and the `Leap_Frog_Position` deferred-sync fix. See
those commits' messages for full narrative detail; this file is the
standing reference for "what's the current performance picture and
what's next," updated as further optimization work lands.

## Benchmark setup

Real system, not a toy topology: `extended_test/ubiquitin_54a8.top` +
`ubiquitin.cnf` — 22712 atoms (775 solute + 7044 SPC solvent). LINCS
(solute, order 4) + M-Shake (solvent) constraints, 2 temperature baths,
isotropic Berendsen pressure coupling, reaction-field electrostatics.
10000 steps, `dt=0.002`, `NTPR`/`NTWX`/`NTWE`/`NSCM` = 1000. Same input
file used for CPU and GPU runs (only the `GPU` block differs). GPU: RTX
50-series class card, confirmed idle (`nvidia-smi`) before each run.

Headline numbers, most recent run:

| | Performance | Wall time (sim) |
|---|---|---|
| CPU (this branch, `USE_CUDA=OFF`) | 2.208 ns/day | 782.7 s |
| CPU (`master`, pre-CUDA-work) | 2.158 ns/day | 800.7 s |
| CPU (this branch, `USE_CUDA=ON`, GPU unused at runtime) | 2.214 ns/day | 780.5 s |
| **GPU (this branch, GPU active)** | **4.88 ns/day** | **343.9 s** |

CPU performance is identical across `master`, CUDA-compiled-but-unused,
and CUDA-disabled builds (within ~2% run-to-run noise) — confirmed by
direct comparison, no CPU-path regression anywhere in this branch's
history. GPU delivers a genuine **~2.2× speedup** over any CPU
configuration. Full TIMING block for the current GPU run (after the
`Leap_Frog_Position` deferred-sync fix, commit `25682dd72`):

```
RemoveCOMMotion                        0.185
Lattice_Shift_Tracker                  0.791
Angle                                 31.345
ImproperDihedral                       7.633
Dihedral                               3.876
Crossdihedral                          0.007
NonBonded                            167.794   (compute forces energies 53.742 / 32.03%, pairlist 113.630 / 67.72%)
MolecularVirial                        4.562
Leap_Frog_Velocity                     9.617
BerendsenThermostat                   29.172
Leap_Frog_Position                     0.082
CUDA_Lincs                             0.313
CUDA_M_Shake                          28.121
TemperatureCalculation                16.768
PressureCalculation                    0.003
BerendsenBarostat                      0.279
```

Previous run (before the `Leap_Frog_Position` fix, commit `48a6decc3`),
for comparison — note `Leap_Frog_Position` 30.863 → 0.082s, and that
~7.7s of that reappeared in `CUDA_M_Shake` (20.372 → 28.121) rather than
vanishing, since `Berendsen_Barostat`'s default `gpu_mirror_touches()`
already forced a publish every step regardless (see "Architecture
direction" below):

```
Leap_Frog_Velocity                     9.591
BerendsenThermostat                   29.059
Leap_Frog_Position                    30.863
CUDA_Lincs                             0.444
CUDA_M_Shake                          20.372
TemperatureCalculation                16.648
PressureCalculation                    0.004
BerendsenBarostat                      0.234
```

Total wall time 369.2s → 343.9s (**~6.9% faster**), muted by the
barostat confound above — expect a larger effect on NVT-only runs or
once a GPU-native barostat exists.

## The core finding: two distinct kinds of GPU cost

Everything in this codebase's GPU port falls into one of two buckets,
and they need different fixes:

1. **Real, embarrassingly-parallel compute** (the nonbonded force
   kernel, bonded-term force kernels, reductions) — genuinely faster on
   GPU, scales with problem size. This is where the actual 2-4.5×
   per-piece speedups come from.
2. **Small/trivial kernels dominated by a fixed CPU↔GPU sync tax.**
   Confirmed empirically: `CUDA_Lincs` has zero blocking calls anywhere
   in its `calculate_interactions()` and costs ~0.04 ms/step; every
   other small algorithm that calls `cudaStreamSynchronize()` (or worse,
   `cudaDeviceSynchronize()`) to read a result back *this step* costs
   roughly a **constant ~2-4 ms/step**, almost totally independent of
   how little data is transferred or how cheap the kernel's own math is.
   This tax is the entire explanation for why `Angle` (28 atoms' worth
   of actual math) costs more wall-clock time per step than the whole
   nonbonded force kernel over ~2.6M pairs.

The fix for bucket 2, proven on `Temperature_Calculation` this session:
**defer the sync** past cheap intervening algorithms until a genuine
consumer needs the result, using the new `Algorithm::finalize_gpu_step()`
/`needs_finalized_gpu_state()` hooks (`algorithm.h`) for `Algorithm`s, or
`Interaction::needs_fresh_cpu_force()` for `Interaction`s (bonded terms/
`NonBonded`). Not every bucket-2 algorithm has a real "wait, is anyone
downstream even going to look at this" gap to exploit — see per-algorithm
notes below.

## Recommendation: what to optimize next

**Top pick — `NonBonded`'s pairlist cost (116.18s, 68% of NonBonded's
own GPU time, the single largest number in the whole profile after the
unavoidable force computation itself).** This is bucket 1, not bucket 2
— real work (candidate search, Morton sort, tile classification over
22712 atoms), already faster than CPU's own pairlist (139.8s), but not
nearly as dramatically improved as the force kernel was (761s → 54s,
4.5× vs pairlist's own 139.8s → 116.2s, only ~1.2×). That gap between
"how much the force kernel improved" and "how much the pairlist
improved" is the biggest remaining opportunity in the codebase, bigger
than the entire bucket-2 cluster combined. It needs real profiling
inside `CUDA_Pairlist_Algorithm_Impl` first (is it rebuild-bound?
classification-bound, which runs every step regardless of the `NSNB`
rebuild cadence? sort-bound?) before proposing a fix — I have not done
that profiling yet, this is a "look here next" flag, not a diagnosed
root cause.

**Lower-risk runner-up, if a quick well-understood win is preferred
instead:** `Berendsen_Thermostat<gpuBackend>`/`NoseHoover_Thermostat<gpuBackend>`
(29.2s) — same bucket-2 shape `Temperature_Calculation` and
`Leap_Frog_Position` were before their fixes, a single-pass reduction
with no sequential-launch-parameter dependency, so the same
`finalize_gpu_step()` deferral pattern should apply directly. Doing
this *before* porting `Berendsen_Barostat` matters: as long as the
barostat forces a full mirror publish every step (see "Architecture
direction" below), any pos/vel-touching bucket-2 fix downstream of it
in the sequence has its win partly absorbed by whichever algorithm
next reads the mirror — `Leap_Frog_Position`'s fix (below) demonstrated
this directly.

(`Leap_Frog_Position`'s own deferred-sync fix — previously listed here
as the runner-up — is now done, commit `25682dd72`. See "Per-algorithm
status" below for its result and the barostat-confound caveat.)

## Per-algorithm status

### GPU-native `Interaction`s (`Forcefield` children)

- **`CUDA_Nonbonded_Interaction`** (`interaction/nonbonded/interaction/`)
  — **solved:** warp-per-tile force kernel rewrite (biggest single win
  of the whole effort, 4.7× on the force-compute sub-line alone), GPU-
  resident force accumulation (writes the mirror directly via
  `mark_gpu_dirty`, no per-call CPU force readback). **Could be solved:**
  pairlist cost, see "Recommendation" above — not yet investigated.
- **`CUDA_Angle_Interaction` / `CUDA_Dihedral_Interaction` /
  `CUDA_Improper_Dihedral_Interaction`** (`interaction/bonded/`) —
  **solved:** GPU-resident force writes (own stream, `atomicAdd` into
  the mirror, `needs_fresh_cpu_force() == false`). **Could be solved:**
  each still does one `cudaStreamSynchronize()` per step to read back
  its own tiny energy/virial buffer — textbook bucket-2 cost (Angle:
  31.5s, ImproperDihedral: 7.1s, Dihedral: 3.9s — roughly proportional
  to term count, consistent with sync tax not compute). Same
  `finalize_gpu_step()`-style deferral as `Temperature_Calculation`
  would apply, deferred until something reads `conf.current().energies`/
  `virial_tensor` (`Molecular_Virial_Interaction` and/or
  `Energy_Calculation`) — not yet done. **`CUDA_Quartic_Bond_Interaction`**
  exists and is ported the same way but wasn't exercised in this
  benchmark (`FORCE` block has `bonds=0`, i.e. bond lengths are
  constrained via LINCS/M-Shake, not force-evaluated) — no data point
  here.
- **`CUDA_Position_Restraint_Interaction`** (`interaction/special/`) —
  **solved:** ported this session, correctness-tested, same GPU-
  resident-force pattern as the other bonded terms. **Not yet
  measured:** this benchmark has no `POSITIONRES` block active, so it
  never runs here — no performance data point yet, only the correctness
  test (`position_restraint_gpu.t.cc`).
- **`Molecular_Virial_Interaction`** — **cannot be solved by porting it
  to GPU** (it's deliberately generic/accelerator-agnostic, correcting
  atomic→molecular virial from the CPU-side force array — that's the
  right design, not a gap). Its real cost (4.3s) is mostly the *force
  publish* it triggers (`needs_fresh_cpu_force() == true` →
  `copy_forces_from_device()`, a genuine one-time-per-step 22712-atom
  device→host element-wise copy). **Could be solved:** that copy loop
  itself is a plain host `for` loop doing per-atom type conversion
  (`configuration_struct.cu`) — replacing it with an async, pinned-
  memory transfer might shrink this, but that's a lower-level I/O
  optimization, not an algorithm redesign, and hasn't been investigated.
- **`Crossdihedral_Interaction`** — **cannot be solved without a real
  GPU port** (no `CUDA_Crossdihedral_Interaction` exists at all — CPU-
  only, unconditionally). Currently irrelevant performance-wise (0.005s,
  this topology has none), but would need a real port (not just a
  deferred-sync fix) if a system with crossdihedral terms mattered.

### `Algorithm`s (`Algorithm_Sequence` members)

- **`Remove_COM_Motion<gpuBackend>`** — **solved:** own stream (not
  device-wide sync), no per-atom host upload loops. **Partially solved:**
  this is a genuinely different shape than the other bucket-2 cases —
  its rotation-removal path needs *two sequential* reduction passes
  where pass 2's kernel launch parameters are pass 1's host-computed
  result (an inherent CPU-in-the-loop dependency, not just "nobody
  downstream needs this yet"). Deferring past other algorithms doesn't
  help the same way it did for `Temperature_Calculation`; a real fix
  would need restructuring the reduction itself (e.g. a single kernel
  that does both passes via a device-side sync primitive, or computing
  the second pass's launch config independent of the first pass's
  result). Not attempted. Currently cheap (0.3s) so low priority.
- **`Lattice_Shift_Tracker`** — **cannot be solved without a real GPU
  port.** Explicitly `static_assert`s against `gpuBackend` in its
  constructor (`lattice_shift.h`) — always runs on CPU regardless of
  accelerator, by original design, not an oversight. Cheap enough
  (0.8s) that porting it likely isn't worth the effort.
- **`Leap_Frog_Velocity<gpuBackend>`** — **solved this session:** fixed
  the accidental coarse-full-resync-every-step bug (`MIRROR_BOX` never
  tracked as fresh, silently forcing a full 8-array upload every call)
  and the `gpu::Configuration` mirror current/old swap gap that was
  hiding behind it. **Could be solved further:** no additional low-
  hanging fruit identified yet at 9.6s; not investigated deeply.
- **`Leap_Frog_Position<gpuBackend>`** — **solved this session:**
  replaced the unconditional `sync_configuration_from_device()` with
  `mark_gpu_dirty(POS|VEL)`; publish now happens lazily via
  `gpu_mirror_touches()`/`flush_gpu_dirty()` for any in-sequence CPU
  consumer, plus a new `io::Out_Configuration::needs_gpu_mirror_flush()`
  gate wired into `program/md.cc` (trajectory writing lives outside
  `Algorithm_Sequence::run()`, so needed its own hook — see
  "Architecture direction" below). Own timer: 30.9s → 0.08s. **Caveat:**
  in *this* benchmark `Berendsen_Barostat` (default `gpu_mirror_touches()
  == MIRROR_ALL`) already forces a full publish every step regardless,
  so ~7.7s of the saved cost reappeared in `CUDA_M_Shake` (20.4s →
  28.1s, the next algorithm that reads position) rather than
  disappearing; net wall time still improved (~6.9%) from the reduced
  transfer volume/timer overhead itself, but the full architectural
  benefit won't show until `Berendsen_Barostat` also stops forcing an
  unconditional publish (NVT-only runs already see it, since there's no
  barostat in the sequence at all).
- **`Berendsen_Thermostat<gpuBackend>` / `NoseHoover_Thermostat<gpuBackend>`**
  (share `thermostat_velocity_scale.h`) — **could be solved:** same
  bucket-2 shape as `Temperature_Calculation` before this session's fix
  — one `cudaStreamSynchronize()` per step to pull back the COM-velocity
  reduction needed for the scale-apply kernel's launch parameters
  (29.1s for Berendsen in this benchmark; NoseHoover shares the same
  code path so likely similar, not separately measured here since this
  benchmark uses Berendsen). Unlike `Remove_COM_Motion`, this reduction
  *is* a single pass (no sequential-launch-parameter dependency), so the
  same `finalize_gpu_step()` deferral pattern should apply directly —
  next-best candidate after `Leap_Frog_Position` if pursuing the
  bucket-2 cluster further.
- **`CUDA_M_Shake`** — **solved:** deferred constraint-error-flags
  (no per-step error-check sync at all, folded into the shared
  end-of-step array). **Could be solved further:** still does one
  `cudaStreamSynchronize()` per step to publish `constraint_force`/
  `virial_tensor` (20.4s) — same category as the bonded terms' energy/
  virial sync, not yet deferred.
- **`CUDA_Lincs`** — **solved, fully.** Zero blocking calls anywhere in
  its `calculate_interactions()` — the reference example of what
  "bucket 2 done right" looks like (0.44s for the same term count class
  as `CUDA_M_Shake`'s constraint work). Nothing left to do here.
- **`Temperature_Calculation<gpuBackend>`** — **solved this session:**
  deferred-sync + merged reduction kernel, 34.5s → 16.6s (2.07×). **Could
  be solved further:** the remaining 16.6s suggests the gap before
  `Energy_Calculation` (its only real consumer) still isn't long enough
  to fully hide the kernel behind other work — `Pressure_Calculation`/
  `Berendsen_Barostat` in between are themselves too cheap to provide
  much cover. Deferring further (if anything later in the sequence could
  serve as the flush point instead) or reducing the kernel count/size
  further are both open, not yet explored.
- **`Pressure_Calculation`** — **solved, trivially:** confirmed it never
  touches the GPU mirror at all (verified by reading the code, not
  guessing) — plain O(9) CPU matrix math from already-CPU-resident
  `virial_tensor`/`kinetic_energy_tensor`. No GPU port needed, ever.
- **`Berendsen_Barostat`** — **could be solved, scope-gated:** confirmed
  free in *this* benchmark (0.23s) because pressure coupling's actual
  per-atom position-scaling loop barely matters at this coupling
  strength/frequency — but it's still a genuine CPU-only, O(num_atoms)
  loop with no GPU backend at all (`Berendsen_Barostat` isn't even
  `Backend`-templated). A real GPU port (isotropic case at minimum) was
  scoped out earlier this session as bigger, multi-mode work not
  justified without a benchmark that actually stresses it — would need
  a system/config with stronger pressure coupling to know if it's worth
  doing.
- **`Energy_Calculation`** — **solved, trivially:** cheap enough it
  doesn't even show as a separate nonzero `TIMING` line. No GPU work
  needed; now doubles as the natural deferred-sync consumer for
  `Temperature_Calculation`'s result.

## Architecture direction: pull-based sync, and a possible step-schedule plan

Working through `Leap_Frog_Position` surfaced a general principle worth
stating explicitly, since it reframes every remaining bucket-2 item
above:

**No GPU-native writer should force a D2H sync. No CPU-only algorithm
should force a GPU-native writer to have synced already.** Instead:
the consumer that actually needs the data pulls it, in whichever
direction it needs (CPU algorithm needing GPU-written data does the
D2H; GPU algorithm needing CPU-written data does the H2D) — never the
producer pushing eagerly "just in case." This is already exactly what
`gpu_mirror_touches()`/`flush_gpu_dirty()`/`configuration_view()` do for
in-sequence algorithms, and what `needs_fresh_cpu_force()` does for
`Interaction`s. `Leap_Frog_Position`'s fix, and the earlier
`io::Out_Configuration::needs_gpu_mirror_flush()` addition, are both
just this same pull-based rule extended to the one consumer that lives
outside `Algorithm_Sequence::run()` (trajectory/checkpoint writing) —
printout doesn't need to join the sequence; it only needs to ask for
the data itself, on its own schedule, which is what it now does.
Read through this lens, `Berendsen_Barostat`'s default
`gpu_mirror_touches() == MIRROR_ALL` is the one remaining *unconditional
push* left in the whole pipeline (every other bucket-2 item pulls
lazily now or is diagnosed-but-not-yet-fixed) — porting it, or at least
narrowing what it touches, is now arguably higher-value than it looked
in isolation, since it's actively muting the payoff of every fix
downstream of it in the sequence (see `Leap_Frog_Position`'s caveat
above).

The other piece raised: **transfer/compute overlap.** Right now every
fix in this file is "avoid the transfer/sync entirely, most steps."
That's the bigger lever, but once bucket-2 syncs are mostly eliminated,
the remaining ones (energy/virial readbacks that gate real downstream
CPU math, e.g. bonded terms' still-per-step sync) are candidates for
overlap rather than elimination — e.g. issuing next step's independent
kernels (velocity update, next pairlist rebuild) on their own streams
before blocking on this step's energy readback, instead of the
current implicit mostly-serial ordering. Not investigated or attempted
yet; flagged here as a second-order optimization after the sync-tax
cluster above is cleared out.

**On a build-time schedule plan:** yes, agreed this is the right
long-term shape, and the metadata to build it already substantially
exists — `gpu_mirror_touches()` per `Algorithm`, `needs_fresh_cpu_force()`/
`is_gpu_native()` per `Interaction`, `needs_finalized_gpu_state()`/
`has_pending_gpu_finalize()` per `Algorithm` — all declarative "what do
I read/write, on which side" facts, currently *consumed* imperatively
(checked fresh at each `apply()` call via `Algorithm_Sequence::run()`'s
loop). A real scheduler would instead *compile* `create_md_sequence.cc`'s
fixed algorithm list once at startup into an explicit per-step plan:
which kernels can issue immediately (no data dependency outstanding),
which syncs are unavoidable and where the best place to put them is
(as late as possible, ideally overlapped with independent GPU work
issued on another stream), and which CPU-only algorithms can run
concurrently with in-flight GPU kernels rather than waiting on the
current mostly-sequential loop. This is a substantially bigger design
than any single fix in this file — a real dependency graph / DAG
scheduler, not just per-algorithm hooks — and hasn't been scoped or
started. Worth a dedicated design pass (probably its own `PLAN.md`
section) once the remaining bucket-2 items are individually fixed and
their sync points are well understood in isolation; premature to build
the general scheduler before that, since the current per-algorithm
fixes are exactly the ground-truth data a scheduler's cost model would
need.

## Verification note

Every "solved" entry above has a passing correctness test
(`ctest`, 30/31 with only the pre-existing `aladip_cuda` perturbation-
gate failure unrelated to any of this) and has been run under
`compute-sanitizer --tool memcheck` with zero errors. "Could be solved"
entries are diagnosed (root cause identified, usually via direct
comparison against `CUDA_Lincs`'s zero-sync reference case) but not yet
implemented.
