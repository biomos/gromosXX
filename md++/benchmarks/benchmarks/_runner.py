"""Run one md++ simulation and return its parsed timing metrics.

Build variants
--------------
Four variants, distinguished by ``GROMOS_VARIANT``. Note that the CPU/GPU
nonbonded kernel is a *runtime* choice made in the imd ``INNERLOOP`` block, not
purely a build-time one, so the ``cuda`` variant is an ordinary OpenMP build
that additionally requests the CUDA kernel per run.

    serial   plain build, one thread
    omp      -DOMP=ON, thread count from the benchmark parameter
    cuda     -DOMP=ON -DCUKERNEL=ON, plus INNERLOOP NTILM=4
    mpi      -DMPI=ON, launched through mpirun; binary is md_mpi

Measurement noise
-----------------
MD benchmarking lives or dies on this, so several things are deliberate:

* trajectory writing is disabled in the imd (see ``_imd``), and the run happens
  in a scratch directory on local disk, so the filesystem is not in the
  measurement;
* the reported metric is ``Wall time simulation``, never ``total`` -- for the
  large system, reading the 29 MB input coordinates takes about four times
  longer than the simulation being measured;
* threads are bound to cores, and the process is pinned to a fixed CPU set, so
  a run is not migrated across sockets mid-measurement;
* each configuration is repeated and the *minimum* is kept. The minimum is the
  right statistic for performance: interference from the rest of the machine
  can only ever make a run slower, so the fastest observation is the closest
  one to the machine's true capability.
"""

import os
import shutil
import subprocess
import sys
import tempfile

from . import _imd, _omd, _systems

VARIANTS = ("serial", "omp", "cuda", "mpi")

# INNERLOOP payloads: NTILM NTILS [NGPUS NDEVG]. NTILM=4 selects the CUDA
# kernel; for CPU variants only two parameters may be present at all.
_INNERLOOP_CPU = [0, 0]
_INNERLOOP_CUDA = [4, 0, 1, 0]


def variant():
    v = os.environ.get("GROMOS_VARIANT", "omp")
    if v not in VARIANTS:
        raise ValueError("unknown GROMOS_VARIANT %r, expected one of %s"
                         % (v, ", ".join(VARIANTS)))
    return v


def repeats():
    return int(os.environ.get("GROMOS_BENCH_REPEATS", "3"))


def binary_name(var):
    return "md_mpi" if var == "mpi" else "md"


def find_binary(var):
    """Locate the md++ binary for this variant.

    asv installs the build into the environment prefix, so inside an asv
    virtualenv ``sys.prefix`` is that environment. ``GROMOS_MD_BINARY`` overrides
    this for standalone use of the harness outside asv.
    """
    name = binary_name(var)

    override = os.environ.get("GROMOS_MD_BINARY")
    if override:
        if not os.path.isfile(override):
            raise RuntimeError("GROMOS_MD_BINARY=%s does not exist" % override)
        if var == "mpi" and os.path.basename(override) != "md_mpi":
            # mpirun happily launches N copies of a non-MPI binary. They do not
            # cooperate: each runs the whole simulation while competing for the
            # same cores, so the result looks like a severe slowdown rather
            # than an error. Refuse rather than record a meaningless number.
            raise RuntimeError(
                "variant 'mpi' needs an MPI build (md_mpi), but "
                "GROMOS_MD_BINARY=%s. Running mpirun on a non-MPI binary "
                "launches independent duplicate simulations and produces "
                "meaningless timings." % override)
        return override

    candidate = os.path.join(sys.prefix, "bin", name)
    if os.path.isfile(candidate):
        return candidate
    found = shutil.which(name)
    if found:
        return found
    raise RuntimeError(
        "cannot find %r (looked in %s and PATH). Set GROMOS_MD_BINARY to run "
        "the harness outside asv." % (name, os.path.dirname(candidate)))


def scratch_root():
    """Local directory for run working files.

    Must be local disk: a network path would put filesystem latency into the
    measurement even with trajectory writing switched off.
    """
    return os.environ.get("GROMOS_BENCH_SCRATCH", "/local/gromos/bench_scratch")


def _cpu_list(n):
    """Fixed CPU set for pinning, honouring any inherited affinity."""
    try:
        available = sorted(os.sched_getaffinity(0))
    except AttributeError:
        available = list(range(os.cpu_count() or 1))
    if n > len(available):
        raise RuntimeError("requested %d threads but only %d CPUs are available"
                           % (n, len(available)))
    return available[:n]


def build_command(var, nthreads, binary, args):
    if var == "mpi":
        return ["mpirun", "-np", str(nthreads), binary] + args
    cmd = [binary] + args
    if os.environ.get("GROMOS_BENCH_PIN", "1") != "0" and shutil.which("taskset"):
        cpus = ",".join(str(c) for c in _cpu_list(nthreads))
        cmd = ["taskset", "-c", cpus] + cmd
    return cmd


def build_env(var, nthreads):
    env = os.environ.copy()
    if var == "mpi":
        # md++ forbids building OMP and MPI together; keep the runtime honest.
        env["OMP_NUM_THREADS"] = "1"
    else:
        env["OMP_NUM_THREADS"] = str(nthreads)
        env["OMP_PROC_BIND"] = "close"
        env["OMP_PLACES"] = "cores"
    return env


def run_once(tier, nthreads, var=None, nstlim=None, keep=False):
    """Run one simulation and return its parsed metrics dict."""
    var = var or variant()
    spec, files = _systems.resolve(tier)
    nstlim = nstlim or spec.get("nstlim_by_variant", {}).get(var, spec["nstlim"])
    binary = find_binary(var)

    os.makedirs(scratch_root(), exist_ok=True)
    workdir = tempfile.mkdtemp(prefix="%s-%s-%dt-" % (tier, var, nthreads),
                               dir=scratch_root())
    try:
        imd_path = os.path.join(workdir, "bench.imd")
        _imd.prepare_benchmark_imd(
            files["imd"], imd_path, nstlim,
            _INNERLOOP_CUDA if var == "cuda" else _INNERLOOP_CPU)

        omd_path = os.path.join(workdir, "bench.omd")
        args = ["@topo", files["top"],
                "@conf", files["cnf"],
                "@input", imd_path,
                "@fin", os.path.join(workdir, "final.cnf")]

        cmd = build_command(var, nthreads, binary, args)
        run_env = build_env(var, nthreads)
        with open(omd_path, "w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                                  env=run_env, text=True)

        try:
            metrics = _omd.parse_file(omd_path)
        except _omd.OmdError as exc:
            stderr = (proc.stderr or "").strip().splitlines()
            tail = "\n  ".join(stderr[-5:]) if stderr else "(no stderr)"
            raise _omd.OmdError(
                "%s\ncommand: %s\nexit code: %d\nstderr tail:\n  %s"
                % (exc, " ".join(cmd), proc.returncode, tail))

        metrics["nstlim"] = float(nstlim)
        metrics["ms_per_step"] = metrics["wall_simulation"] / nstlim * 1000.0
        _trace(tier, var, nthreads, nstlim, metrics, run_env)
        return metrics
    finally:
        if not keep:
            shutil.rmtree(workdir, ignore_errors=True)


def _trace(tier, var, nthreads, nstlim, metrics, run_env=None):
    """Append one line per simulation to $GROMOS_BENCH_TRACE, if set.

    Useful for confirming how many times a configuration actually ran, and
    under what conditions -- a benchmark that silently runs more often than
    expected, or on fewer cores than requested, is hard to spot from the
    reported numbers alone.
    """
    path = os.environ.get("GROMOS_BENCH_TRACE")
    if not path:
        return
    try:
        affinity = len(os.sched_getaffinity(0))
    except AttributeError:
        affinity = -1
    with open(path, "a") as fh:
        fh.write("pid=%d %s var=%s threads=%s omp=%s affinity=%d nstlim=%d "
                 "sim=%.3f ms/step=%.3f\n"
                 % (os.getpid(), tier, var, nthreads,
                    (run_env or {}).get("OMP_NUM_THREADS"), affinity, nstlim,
                    metrics["wall_simulation"], metrics["ms_per_step"]))


def run_best(tier, nthreads, var=None, nstlim=None, n=None):
    """Repeat a configuration and keep the best (fastest) observation.

    All metrics are taken from the single repeat with the lowest simulation
    wall time, rather than reducing each metric independently -- that keeps the
    per-algorithm breakdown internally consistent and summing to the total.
    """
    n = repeats() if n is None else n
    best = None
    for _ in range(max(1, n)):
        metrics = run_once(tier, nthreads, var=var, nstlim=nstlim)
        if best is None or metrics["wall_simulation"] < best["wall_simulation"]:
            best = metrics
    return best
