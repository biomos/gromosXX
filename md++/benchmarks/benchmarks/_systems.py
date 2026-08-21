"""Benchmark system registry.

Three equilibrated protein-in-solvent systems spanning roughly an order of
magnitude in size. All three use the same reaction-field setup (0.8/1.4 nm
cutoffs, RCRF eps 62, SHAKE with NTC=3, pairlist rebuilt every 5 steps), so a
comparison across tiers isolates system size rather than method.

The systems originate from a curated protein set that lives on shared storage,
but the benchmark harness never reads that storage: ``tools/stage_systems.py``
copies what is needed into a local cache once, and everything here resolves
against that cache. See ``README.md`` for the rationale.
"""

import json
import os

# Where the staged systems live. Must be local storage on the benchmark machine:
# a network mount would make the timeline measure the network.
DEFAULT_CACHE = "/local/gromos/bench_systems"

# Upstream source, read only by tools/stage_systems.py.
DEFAULT_SOURCE = "/pool/discus_2c/prot_set_mm/gromos/RF"

MANIFEST_NAME = "manifest.json"

# NSTLIM is tuned so a run takes roughly 30 s at 4 OpenMP threads. Cost scales
# close to linearly with atom count (measured: 8.1 / 33.9 / 136.8 ms per step
# for tiny / medium / large). `tools/calibrate.py` re-measures and suggests
# updated values.
#
# A single NSTLIM per system keeps a tier directly comparable across build
# variants, at the cost of uneven durations: serial runs take ~2.75x longer than
# 4-thread ones and CUDA runs about half as long. If that becomes a problem, a
# system may carry an optional per-variant override::
#
#     "nstlim_by_variant": {"serial": 1500, "cuda": 8000}
#
# which _runner consults before falling back to "nstlim". Note that changing
# either value breaks comparability with results already recorded.
SYSTEMS = {
    "tiny": {
        "pdb": "2ovn",
        "atoms": 11904,
        "nstlim": 4000,
    },
    "medium": {
        "pdb": "1aki",  # hen egg-white lysozyme
        "atoms": 45530,
        "nstlim": 1000,
    },
    "large": {
        "pdb": "1tua",
        "atoms": 172541,
        "nstlim": 250,
    },
}

TIERS = ("tiny", "medium", "large")


def cache_root():
    """Directory holding the staged systems."""
    return os.environ.get("GROMOS_BENCH_SYSTEMS", DEFAULT_CACHE)


def source_root():
    """Upstream directory to stage from. Only used by tools/stage_systems.py."""
    return os.environ.get("GROMOS_BENCH_SYSTEMS_SRC", DEFAULT_SOURCE)


def source_files(pdb):
    """Files to copy out of the upstream set for one system.

    ``prep_sys_fin.cnf`` is the equilibrated coordinate file. It is the only
    starting structure present for every system in the upstream set -- the
    per-nanosecond production snapshots are incomplete (1aki, for instance, has
    no md_99.cnf), so relying on them would silently break some systems.
    """
    return {
        "top": os.path.join(source_root(), pdb, "prep_sys_0", "%s_ions.top" % pdb),
        "cnf": os.path.join(source_root(), pdb, "prep_sys_0", "prep_sys_fin.cnf"),
        "imd": os.path.join(source_root(), pdb, "prod_md_0", "md_100.imd"),
    }


def staged_files(pdb):
    """Paths to one system's files inside the local cache."""
    base = os.path.join(cache_root(), pdb)
    return {
        "top": os.path.join(base, "%s_ions.top" % pdb),
        "cnf": os.path.join(base, "prep_sys_fin.cnf"),
        "imd": os.path.join(base, "md_100.imd"),
    }


def resolve(tier):
    """Return ``(spec, files)`` for a tier, or raise if it is not staged."""
    try:
        spec = SYSTEMS[tier]
    except KeyError:
        raise KeyError("unknown system tier %r, expected one of %s"
                       % (tier, ", ".join(TIERS)))
    files = staged_files(spec["pdb"])
    missing = [p for p in files.values() if not os.path.isfile(p)]
    if missing:
        raise RuntimeError(
            "system %r (%s) is not staged -- missing:\n  %s\n"
            "Run: python tools/stage_systems.py"
            % (tier, spec["pdb"], "\n  ".join(missing)))
    return spec, files


def load_manifest():
    path = os.path.join(cache_root(), MANIFEST_NAME)
    if not os.path.isfile(path):
        raise RuntimeError(
            "no %s in %s -- run: python tools/stage_systems.py"
            % (MANIFEST_NAME, cache_root()))
    with open(path) as fh:
        return json.load(fh)
