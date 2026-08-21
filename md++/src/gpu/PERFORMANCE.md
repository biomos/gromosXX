# PERFORMANCE.md — CUDA port performance state

Snapshot as of the latest commit (branch `cuda_claude_rf_excluded`).
Written after a benchmarking/optimization pass covering the async
stream-ordered redesign, the `POSITIONRES`-dropped-force bug fix, the
mixed-CPU/GPU warning mechanism, the `Temperature_Calculation` and
`Leap_Frog_Position` deferred-sync fixes, `Berendsen_Barostat`'s GPU
port plus `RemoveCOMMotion`/`Lattice_Shift_Tracker`'s deferred-sync
fixes, and — this pass — switching force/constraint-force/virial/
kinetic-energy/pressure-tensor accumulation to double precision
(`FPH`), fixing `Berendsen_Thermostat`/`NoseHoover_Thermostat`'s
per-step sync tax, and diagnosing (not yet fixing) the root cause of
the "residual resync tug-of-war" from the previous pass — see
"Diagnosed: Unified Memory is the real cost" below. See individual
commits' messages for full narrative detail; this file is the standing
reference for "what's the current performance picture and what's
next," updated as further optimization work lands.

## Benchmark setup

Real system, not a toy topology: `extended_test/ubiquitin_54a8.top` +
`ubiquitin.cnf` — 22712 atoms (775 solute + 7044 SPC solvent). LINCS
(solute, order 4) + M-Shake (solvent) constraints, 2 temperature baths,
isotropic Berendsen pressure coupling, reaction-field electrostatics.
10000 steps, `dt=0.002`, `NTPR`/`NTWX`/`NTWE`/`NSCM` = 1000. Same input
file used for CPU and GPU runs (only the `GPU` block differs). GPU: RTX
50-series class card, confirmed idle (`nvidia-smi`) before each run.

Headline numbers, most recent run (GPU with `FPH`-precision force/
virial accumulation and the `Berendsen_Thermostat` sync fix; CPU
run on the same binary, `NTGPU=0`, for a true apples-to-apples
comparison):

| | Performance | Wall time (sim) |
|---|---|---|
| CPU (this branch, `USE_CUDA=OFF`, or `NTGPU=0` on a `USE_CUDA=ON` build) | ~2.2-2.25 ns/day | ~769-783 s |
| **GPU (this branch, GPU active, `FPH` accumulation)** | **4.60 ns/day** | **375.7 s** |

GPU delivers a genuine **~2.0-2.1× speedup** over CPU. This is down
from the ~2.2-2.3× measured before this pass — see "The `FPH`
precision tradeoff" below for why, and why it was still the right
call. Full TIMING block for the current GPU run:

```
RemoveCOMMotion                        0.158
Lattice_Shift_Tracker                  0.098
Angle                                  9.116
ImproperDihedral                       6.403
Dihedral                               3.898
Crossdihedral                          0.003
NonBonded                            197.343   (compute forces energies 73.521 / 37.26%, pairlist 123.419 / 62.54%)
MolecularVirial                        4.486
Leap_Frog_Velocity                    24.380
BerendsenThermostat                    0.285
Leap_Frog_Position                     0.058
CUDA_Lincs                             0.275
CUDA_M_Shake                          36.651
TemperatureCalculation                31.933
PressureCalculation                    0.001
BerendsenBarostat                      0.068
```

Same run, CPU-only (`NTGPU=0`, same binary/topology/config), for direct
comparison — note the completely different cost shape (`NonBonded` is
97% of CPU time; every "expensive on GPU" bucket-2 line above costs
essentially nothing on CPU, since CPU has no fixed per-call sync tax to
pay):

```
RemoveCOMMotion                        0.002
Lattice_Shift_Tracker                  0.750
Angle                                  0.358
ImproperDihedral                       0.246
Dihedral                               0.635
Crossdihedral                          0.002
NonBonded                            748.880   (shortrange solvent-solvent 313.763 / 41.90%, longrange solvent-solvent 263.991 / 35.25%, pairlist 137.355 / 18.34%)
MolecularVirial                        3.939
Leap_Frog_Velocity                     0.511
BerendsenThermostat                    0.220
Leap_Frog_Position                     0.362
Lincs                                  1.325   (solute 1.171 / 88.41%)
M_Shake                                7.442   (solvent 6.782 / 91.14%)
TemperatureCalculation                 1.751
PressureCalculation                    0.002
BerendsenBarostat                      0.174
```

`NonBonded` is where GPU wins decisively: 748.9s → 197.3s, a genuine
**3.8× speedup** on real embarrassingly-parallel compute. Every other
line is the "two buckets" model in miniature: cheap on CPU (no sync tax
to pay), 10-100× more expensive in *absolute* terms on GPU purely from
CPU↔GPU synchronization overhead — individually a rounding error
against `NonBonded`'s 197s, but collectively enough to dilute a 3.8×
per-kernel win down to ~2.0× overall.

## The `FPH` precision tradeoff

`force`, `constraint_force`, and the mirror's `virial_tensor`/
`kinetic_energy_tensor`/`pressure_tensor` were `FPL` (low precision,
float under this build's `FP_PRECISION=2`) — including the shared
accumulator every bonded/nonbonded/constraint kernel `atomicAdd`s into.
Summing hundreds of small float contributions per atom (nonbonded
pairs, bonded terms, constraint corrections) into one shared total is
exactly where float accumulation error compounds. Switched every such
accumulator to `FPH` (double under `FP_PRECISION` 2/3): per-pair/
per-term *compute* stays `FPL` (that's the bulk of the FLOP volume, and
individual-term precision loss is fine), only the `atomicAdd`
*destination* widened, with an explicit cast at the call site — see
`gpu::CuVArrayH` (`cuvector.h`) and every bonded/nonbonded/constraint
kernel's `force`/`constraint_force` parameter type.

Cost: real, and it's exactly what you'd expect — `NonBonded`'s
"compute forces energies" sub-line went 53.9s → 80.8s(ish, run-to-run)
before the thermostat fix rebalanced other lines; double-precision
`atomicAdd` runs at roughly half the throughput of float on this GPU
generation (fewer FP64 ALUs than FP32), and `NonBonded` is by far the
highest-volume accumulator target in the system (millions of pairwise
contributions converging onto 22712 atoms' force buffer every step).
Net effect on overall GPU-vs-CPU speedup: ~2.2-2.3× → ~2.0-2.1×. Judged
worth it: this is a correctness fix (eliminating float-accumulation
error in the physically-important quantities: total force, virial,
kinetic energy, pressure), not a style preference, and GPU is still a
clear win over CPU even with the cost paid.

## Diagnosed: Unified Memory is the real cost behind bucket 2

The previous pass's "residual resync tug-of-war" (whichever algorithm
did the first POS resync of a step paid ~15s/10000-steps, and fixing
one algorithm just moved the cost to whichever became "first" next)
had two open hypotheses. It's now understood: **`math::CuVArray`/
`gpu::cuvector` (`cuvector.h`) are backed by `cudaMallocManaged` --
CUDA Unified/Managed Memory -- not plain `cudaMalloc`.** Every
`copy_pos_vel_to_device()`/`copy_to_device()`/`copy_pos_vel_from_
device()`-style routine (`configuration_struct.cu`) is a *host*
`for`-loop touching these arrays element-by-element. Unified Memory
pages migrate on demand based on which processor last touched them --
if a page was last written by a GPU kernel and a host loop then reads/
writes it, the CUDA driver faults, waits for outstanding GPU work on
that memory to finish, and physically migrates the page over PCIe/
NVLink. That's a real, host-blocking stall, entirely invisible in the
source (no explicit `cudaMemcpy`/`cudaDeviceSynchronize` call to point
at) and **not tied to which algorithm "logically" needed the sync** --
only to whichever host loop happens to touch a GPU-resident page
first. This is exactly why the cost moved between `Lattice_Shift_
Tracker` and `Leap_Frog_Velocity` depending on `Energy_Calculation`'s
mask in the previous pass: the mask changed *which* host loop was
first to touch stale-on-GPU pages, not whether a migration had to
happen.

**The real fix, not yet implemented:** route the coarse `copy_*`
paths through a persistent pinned (page-locked, `cuhvector`/
`CuHAllocator`, already exists in this codebase) staging buffer instead
of touching the Unified Memory arrays directly from the host. Splits
into two independently-cheap steps: GPU array <-> pinned buffer (same
precision both sides, so a pure `cudaMemcpyAsync` -- real DMA at
near-peak bandwidth, no page faults, genuinely async) and pinned buffer
<-> CPU `double` array (pure host-to-host cast, no device involvement,
no driver synchronization, just a cache-friendly loop). Also worth
reconsidering: do `pos`/`vel`/`force` need to be Unified Memory at all,
if nothing legitimately needs to read/write them directly from the
host anymore once the pinned-staging path exists? Not scoped or
started this session -- flagged here as the next concrete piece of
work, now that the "why" is understood rather than guessed at.

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

**Second pick, if the pairlist investigation isn't the priority:** the
pinned-staging transfer fix ("Diagnosed: Unified Memory is the real
cost" above) — this is the one fix that would validate or invalidate a
lot of this session's "moves, doesn't vanish" findings at once
(`Leap_Frog_Velocity`, `CUDA_M_Shake`, `Temperature_Calculation` all
show the same shape). Higher confidence than the pairlist item since
the root cause is now understood, not just flagged.

(`Leap_Frog_Position`, `Berendsen_Barostat`/`RemoveCOMMotion`/`Lattice_
Shift_Tracker`, and `Berendsen_Thermostat`/`NoseHoover_Thermostat`'s
deferred-sync fixes are all done. See "Per-algorithm status" below.)

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

- **`Remove_COM_Motion<gpuBackend>`** — **solved this session (sync
  deferral):** own stream (not device-wide sync), no per-atom host
  upload loops (both pre-existing), plus `gpu_mirror_touches()` narrowed
  to 0 and `sync_configuration_from_device()` replaced with
  `mark_gpu_dirty()` — same pattern as `Leap_Frog_Position`. Matters
  disproportionately here: this is the *first* algorithm in
  `create_md_sequence.cc`'s sequence, so its old default `MIRROR_ALL`
  forced a full mirror round-trip at the very top of every step,
  regardless of `comtransrot`'s skip-step cadence. Now effectively free
  (0.25s, all in the ~1% of steps that do real removal work).
  **Partially solved still:** the rotation-removal path's inherent two-
  sequential-reduction-pass shape (pass 2's launch parameters are pass
  1's host-computed result) is unchanged — not a sync-deferral problem,
  a genuine algorithmic dependency; only matters on real-removal steps,
  low priority given those are now rare and each one is cheap anyway.
- **`Lattice_Shift_Tracker<gpuBackend>`** — **solved this session (real
  GPU port):** one thread per chargegroup, reuses `gpu::Topology`'s
  existing chargegroup-offset array (`gpu::Periodicity<B>::put_into_box()`
  for the wrap, self-contained kernel for the shift bookkeeping — see
  `lattice_shift_kernels.cu`). `gpu_mirror_touches()==0`, deferred
  publish, same as every other GPU-native writer. Own timer: 0.79s
  (CPU) → 0.10s. Needed two follow-up fixes to actually land clean, both
  documented in "Architecture direction" below: (1) the new
  `MIRROR_LATTICE_SHIFT` field was initially included in `MIRROR_ALL`,
  which made every ordinary algorithm's default touch force a real
  re-upload (0.8s → 17.5s, worse than the CPU version, before this was
  caught and fixed); (2) once excluded from `MIRROR_ALL`, the field had
  no "starts fresh" point at all, so its producer-events list grew by
  one entry every step forever — an O(steps²) cost invisible under
  ~2000 steps, dominant by 10000 (11.5-17.5s) — fixed by adding
  `CudaManager::clear_stale_producer_events()`. Only vacuum/rectangular
  boundary supported (matches `CUDA_Pairlist_Algorithm`'s own scope);
  `init()` hard-errors on triclinic/truncoct rather than silently
  wrapping wrong. `full_anisotropic`-equivalent generality wasn't the
  issue here (chargegroup wrap doesn't depend on pressure-coupling
  mode) — this is purely a boundary-condition scope limit.
- **`Leap_Frog_Velocity<gpuBackend>`** — (earlier session) fixed the
  accidental coarse-full-resync-every-step bug (`MIRROR_BOX` never
  tracked as fresh) and the mirror current/old swap gap hiding behind
  it; stable at ~9.6s for several benchmark runs after that.
  **Regressed, not yet fixed:** 9.6s → ~24s once `Energy_Calculation`
  was narrowed (a later pass). Root cause of the *magnitude* is the
  Unified Memory page-migration cost (see "Diagnosed" section above);
  a targeted attempt to fix the *trigger* (narrowing `Forcefield`'s
  mask, since it's the algorithm right before this one that still
  forces a POS/VEL invalidation) broke `CUDA_Shake` correctness and was
  reverted — see "`Leap_Frog_Velocity`/`Forcefield`" below. Two
  independent pieces both need to land before this is actually fixed:
  the pinned-staging transfer fix (removes the *cost* of whichever
  resync happens) and understanding the `CUDA_Shake` dependency (removes
  the *unnecessary* resync itself).
- **`Leap_Frog_Position<gpuBackend>`** — **solved:** replaced the
  unconditional `sync_configuration_from_device()` with
  `mark_gpu_dirty(POS|VEL)`; publish now happens lazily via
  `gpu_mirror_touches()`/`flush_gpu_dirty()` for any in-sequence CPU
  consumer, plus `io::Out_Configuration::needs_gpu_mirror_flush()`
  gating trajectory writes in `program/md.cc` (outside `Algorithm_
  Sequence::run()`, needed its own hook). Own timer: 30.9s → 0.08s,
  stable across this session's later changes too.
- **`Berendsen_Thermostat<gpuBackend>` / `NoseHoover_Thermostat<gpuBackend>`**
  (share `thermostat_velocity_scale.h`) — **solved this session:** the
  `cudaStreamSynchronize()` that used to pull the COM-velocity reduction
  back to the host (just to compute `com_v_per_group[g] = sums/mass`, a
  trivial O(num_groups) divide) is gone — that divide now happens in a
  new tiny kernel (`gpu::launch_group_com_velocity`, `temperature_
  kernels.cu`), so reduction → divide → scale-apply is three kernels on
  one stream with zero host round trips, same "merge the trivial
  postprocess into a kernel" pattern already proven on `Temperature_
  Calculation`. Own timer: 28.6s → 0.29s (**~100×**, now `CUDA_Lincs`-
  level). NoseHoover shares the identical code path, so should see the
  same fix, not separately re-measured.
- **`CUDA_M_Shake`** — **solved:** deferred constraint-error-flags
  (no per-step error-check sync at all, folded into the shared
  end-of-step array). **Could be solved further:** still does one
  `cudaStreamSynchronize()` per step to publish `constraint_force`/
  `virial_tensor` (currently ~27s, drifted up from 20.4s across this
  session's changes along with the same residual-resync pattern noted
  under `Leap_Frog_Velocity` — not separately diagnosed) — same
  category as the bonded terms' energy/virial sync, not yet deferred.
- **`CUDA_Lincs`** — **solved, fully.** Zero blocking calls anywhere in
  its `calculate_interactions()` — the reference example of what
  "bucket 2 done right" looks like (0.44s for the same term count class
  as `CUDA_M_Shake`'s constraint work). Nothing left to do here.
- **`Temperature_Calculation<gpuBackend>`** — deferred-sync + merged
  reduction kernel from an earlier session (34.5s → 16.6s). **Regressed
  this session, same root cause as `Leap_Frog_Velocity`:** 16.6s → 32.8s
  after `Energy_Calculation`'s narrowing — this is its own deferred-
  finalize consumer, so the two are directly coupled; not separately
  bisected from the `Leap_Frog_Velocity` case, presumed same underlying
  issue. See "residual resync tug-of-war" below.
- **`Pressure_Calculation`** — **solved this session (narrowed):**
  confirmed it only reads `conf.old().virial_tensor`/
  `kinetic_energy_tensor` and writes `pressure_tensor` — plain CPU
  matrices, never GPU-mirror pos/vel/force/box. `gpu_mirror_touches()`
  narrowed from the base class default `MIRROR_ALL` to `MIRROR_VIRIAL`
  (the one field it genuinely needs flushed, since GPU constraint
  algorithms write constraint-virial there). This is the change that
  first surfaced the "any default-`MIRROR_ALL` algorithm forces the
  round-trip regardless of what upstream deferred" pattern — see
  "Architecture direction" below.
- **`Berendsen_Barostat<gpuBackend>`** — **solved this session (real GPU
  port):** isotropic/anisotropic/semi-anisotropic all reduce to "scale
  every position by a fixed 3x3 matrix," one generic kernel covers all
  three (`berendsen_barostat_kernels.cu`); mu/box computed on the host
  (O(9), cheap, identical to the CPU formulas). `full_anisotropic`
  refused with an explicit error rather than silently replicating a
  pre-existing CPU-only bug (the CPU reference's `case pcouple_full_
  anisotropic` block has no `break;` before `case pcouple_semi_
  anisotropic`, so it silently re-applies semi-anisotropic scaling on
  top of its own result — not this branch's call to fix as a side
  effect of a performance patch). Own timer: 0.28s → 0.08s.
- **`Energy_Calculation`** — **narrowed this session, net-positive but
  not clean:** confirmed it only touches `conf.old().energies`
  (already CPU-resident) and `conf.current().averages` — narrowed
  `gpu_mirror_touches()` from `MIRROR_ALL` to 0. Fixed `Lattice_Shift_
  Tracker`'s cost (see above) but shifted a residual POS-resync cost
  onto `Leap_Frog_Velocity`/`Temperature_Calculation` instead of
  eliminating it — bisected and measured (see "residual resync
  tug-of-war" below): keeping the narrowing gives 345.3s total,
  reverting it gives 360.0s, so it's the better of the two known
  options, but neither is a clean win. Still doubles as the deferred-
  sync consumer for `Temperature_Calculation`'s result.

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

### Two case studies from this pass: what "narrow the mask" can still get wrong

**`MIRROR_LATTICE_SHIFT` in `MIRROR_ALL` (found and fixed).** Adding a
new mirror-tracked field and including it in `MIRROR_ALL` "to be safe"
(matching `MIRROR_CONSTRAINT_FORCE`/`MIRROR_VIRIAL`'s precedent) is
*not* automatically safe the way it looks: `MIRROR_ALL` membership
means every ordinary algorithm's default touch **invalidates** the
field, and invalidation isn't a free no-op the way a same-value flush
is — the next read has to pay a real re-upload. A field genuinely
read by other code (like `MIRROR_VIRIAL`, which `Pressure_Calculation`
needs) belongs in `MIRROR_ALL`; a field with exactly one writer and no
other reader anywhere in the codebase (like `MIRROR_LATTICE_SHIFT`)
does not, no matter how "safe" excluding it looks on first read of the
diff. Lesson: check for *any* other consumer, not just "does this look
like the same shape as an existing bit," before deciding `MIRROR_ALL`
membership.

**Producer-events list with no "starts fresh" point (found and
fixed).** Every mirror-tracked field needs some point in its lifecycle
where `field_producer_events[bit]` gets cleared — normally either
`invalidate_gpu_mirror()` (generic, via `MIRROR_ALL` default touches)
or an explicit per-algorithm equivalent (`zero_mirror_force()` for
`MIRROR_FORCE`). A field excluded from `MIRROR_ALL` specifically
*because* nothing else should invalidate it (the `MIRROR_LATTICE_SHIFT`
fix above) loses that generic clearing point as a side effect — the
writer must supply its own, or the list grows by one event every
step forever, and every future `cudaStreamWaitEvent()` pass over that
list gets slower every step (O(steps²) total). Added `CudaManager::
clear_stale_producer_events()` (clears events without touching
freshness/dirty bits, unlike `invalidate_gpu_mirror()`) for exactly
this case; any future genuinely-single-writer field excluded from
`MIRROR_ALL` will need the same treatment.

### The residual resync tug-of-war — root cause now diagnosed

See "Diagnosed: Unified Memory is the real cost behind bucket 2" above
for the full explanation — `cudaMallocManaged`-backed arrays touched by
host `for`-loops fault-and-migrate on demand, invisibly and
independent of which algorithm "logically" triggered the resync. This
is why fixing one algorithm's sync tax kept relocating the cost to
whichever became "first to touch a stale page" next
(`Lattice_Shift_Tracker` ↔ `Leap_Frog_Velocity`/`Temperature_
Calculation`, depending on `Energy_Calculation`'s mask) rather than
eliminating it. The fix (pinned staging buffers, not yet implemented)
is scoped above; once it lands, re-measure whether `Energy_
Calculation`'s narrowing is still the right call, or whether the
tug-of-war disappears entirely and the mask choice stops mattering.

### `Leap_Frog_Velocity`/`Forcefield` — narrowing attempted, reverted, open

`Forcefield::gpu_mirror_touches()` (`forcefield.h`) still returns
`MIRROR_ALL & ~(FORCE|VIRIAL)`, i.e. `POS|VEL|BOX|CONSTRAINT_FORCE` —
despite its own doc comment stating outright that Forcefield never
writes any of those. Tried narrowing to `0u` entirely (it's the very
next default-`MIRROR_ALL` algorithm to touch POS/VEL before `Leap_Frog_
Velocity`, which needs POS+VEL+FORCE immediately after — a plausible
source of exactly the "unnecessary resync" pattern this whole session
has been chasing). **Broke `CUDA_Shake` correctness**: `ubiquitin_gpu`
failed with "SHAKE error, vectors orthogonal" by step 3. Root cause not
understood: `CUDA_Shake::apply()` (`cuda_shake.cc`) uploads `conf.
current()/old().pos` directly from the CPU array (`vec3_upload_fpl`,
bypassing `configuration_view()` entirely) and has its own default
`gpu_mirror_touches() == MIRROR_ALL`, which *should* independently
guarantee that CPU array is fresh before it reads it, regardless of
`Forcefield`'s mask — but bisection is unambiguous (bisected against
the companion `Berendsen_Thermostat` fix too, to rule out an
interaction effect: `Forcefield`'s narrowing alone, with the thermostat
fix reverted, still reproduces the failure). Reverted rather than
shipped with an unexplained correctness dependency. **Next step**: trace
exactly which mirror-freshness state differs between the two mask
values at the point `CUDA_Shake` reads `conf.current().pos` — likely
needs instrumenting `gpu_fresh_fields`/`gpu_dirty_fields` directly
rather than reasoning about it, since the reasoning above says this
*shouldn't* depend on `Forcefield` at all and evidently does.

## Verification note

Every "solved" entry above has a passing correctness test
(`ctest`, 31/32 with only the pre-existing `aladip_cuda` perturbation-
gate failure unrelated to any of this) and has been run under
`compute-sanitizer --tool memcheck` with zero errors — including the 13
force-touching tests re-verified after the `FPH` precision switch
(every bonded/nonbonded/constraint GPU test plus `ubiquitin_gpu`).
"Could be solved" entries are diagnosed (root cause identified, usually
via direct comparison against `CUDA_Lincs`'s zero-sync reference case)
but not yet implemented. Step-0 `E_Total` (`-2.4963e+05`) has been
checked unchanged across every benchmark run referenced in this file,
CPU and GPU alike (confirmed identical in the direct CPU-vs-GPU
comparison run this pass).

The `Forcefield`/`Leap_Frog_Velocity` narrowing attempt (this pass) is
the one exception to "every change here shipped tested": it was
implemented, found to break `CUDA_Shake` on `ubiquitin_gpu`, and
reverted before commit — the code as committed does not include that
narrowing. See "`Leap_Frog_Velocity`/`Forcefield`" above.
