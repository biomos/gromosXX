"""Measure per-step cost and suggest NSTLIM values.

NSTLIM has to be chosen per system so that a run is long enough to average out
noise but short enough to be worth running per commit. This measures the real
cost with a short probe run and reports the NSTLIM that would hit a target
duration.

Doubles as the standalone smoke test for the harness: it exercises system
resolution, imd rewriting, launching and omd parsing without involving asv.

    python tools/calibrate.py --system tiny --threads 4
    GROMOS_MD_BINARY=../BUILD/program/md python tools/calibrate.py --all
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from benchmarks import _runner, _systems


def check_governor():
    """Warn if the CPU is free to change frequency underneath the benchmark."""
    path = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"
    try:
        with open(path) as fh:
            gov = fh.read().strip()
    except OSError:
        return
    if gov != "performance":
        print("warning: CPU governor is %r, not 'performance'. Frequency "
              "scaling adds run-to-run noise and can mask real regressions.\n"
              "         sudo cpupower frequency-set -g performance\n" % gov,
              file=sys.stderr)


def calibrate(tier, threads, probe, target, var):
    spec = _systems.SYSTEMS[tier]
    print("%-6s %-5s %7d atoms  variant=%s threads=%d"
          % (tier, spec["pdb"], spec["atoms"], var, threads), flush=True)

    m = _runner.run_once(tier, threads, var=var, nstlim=probe)
    ms = m["ms_per_step"]
    suggested = max(10, int(round(target * 1000.0 / ms / 10.0)) * 10)

    print("  probe            %d steps" % probe)
    print("  init             %8.3f s" % m["wall_initialisation"])
    print("  simulation       %8.3f s" % m["wall_simulation"])
    print("  per step         %8.3f ms" % ms)
    print("  ns/day           %8.3f" % m["ns_per_day"])
    print("  configured       NSTLIM=%d  (~%.1f s)"
          % (spec["nstlim"], spec["nstlim"] * ms / 1000.0))
    print("  suggested        NSTLIM=%d  (~%.0f s target)" % (suggested, target))

    top = sorted(((v, k) for k, v in m.items() if k.startswith("timer.")),
                 reverse=True)[:5]
    print("  dominant timers: "
          + ", ".join("%s %.2fs" % (k.split(".", 1)[1], v) for v, k in top))
    print(flush=True)
    return suggested


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", default=None, help="tier to calibrate")
    ap.add_argument("--all", action="store_true", help="calibrate every tier")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--variant", default=None,
                    help="serial|omp|cuda|mpi (default: $GROMOS_VARIANT or omp)")
    ap.add_argument("--probe", type=int, default=100,
                    help="steps for the probe run (default: 100)")
    ap.add_argument("--target", type=float, default=30.0,
                    help="target run duration in seconds (default: 30)")
    args = ap.parse_args()

    if not args.all and not args.system:
        ap.error("pass --system TIER or --all")

    check_governor()
    var = args.variant or _runner.variant()
    tiers = list(_systems.TIERS) if args.all else [args.system]

    print("binary: %s\n" % _runner.find_binary(var), flush=True)
    for tier in tiers:
        calibrate(tier, args.threads, args.probe, args.target, var)


if __name__ == "__main__":
    main()
