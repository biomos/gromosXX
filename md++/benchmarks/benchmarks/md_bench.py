"""Plain MD benchmarks for md++.

Each tracked quantity is a ``track_*`` function rather than an asv ``time_*``
benchmark. asv's timing machinery is built for microbenchmarks: it would re-run
each configuration many times and measure the whole process, so a multi-minute
MD run would be timed together with process startup, topology parsing and the
input coordinate read. md++ already reports its own wall times and a
per-algorithm breakdown, so the far cheaper and more informative approach is to
run once and scrape those numbers.

That is what ``setup_cache`` is for: asv calls it once per benchmark class and
hands its return value to every benchmark method. All the simulations happen
there, and each ``track_*`` is then a dictionary lookup. Doing the run in
``setup`` instead would re-run the simulation once per metric.

Note that ns/day is deliberately *not* the headline metric. asv's regression
detection assumes lower is better, so on a higher-is-better series every genuine
speedup would be reported as a regression. Seconds are tracked as the primary
signal; ns/day is carried alongside for readability.
"""

import os

from . import _runner, _systems

# Thread/rank counts. Kept to powers of two up to 16 on the assumption of a
# machine with at least that many physical cores; _runner refuses counts that
# exceed the available CPU affinity rather than silently oversubscribing.
THREADS = [1, 2, 4, 8, 16]

TIERS = list(_systems.TIERS)


def _selected_threads():
    """Thread counts to actually run, which may be a subset of the grid."""
    override = os.environ.get("GROMOS_BENCH_THREADS")
    if override:
        return [int(t) for t in override.replace(",", " ").split()]
    return list(THREADS)


def _selected_tiers():
    """Systems to actually run, which may be a subset of the grid."""
    override = os.environ.get("GROMOS_BENCH_TIERS")
    if override:
        return [t.strip() for t in override.replace(",", " ").split()]
    return list(TIERS)


def _supported(var, nthreads, available):
    """Whether a (variant, thread count) pair is worth running.

    The serial build has no thread parallelism at all, so anything above one
    thread would silently repeat the same measurement.
    """
    if nthreads > available:
        return False
    if var == "serial" and nthreads != 1:
        return False
    return True


class PlainMD(object):
    """Plain MD on protein-in-solvent systems of increasing size."""

    # The parameter grid is fixed, never derived from the environment. asv keys
    # results by position in this grid and records it in benchmarks.json, so a
    # grid that changed shape between invocations would silently orphan every
    # result recorded under a different one -- `asv publish` then emits empty
    # graphs while reporting success.
    #
    # GROMOS_BENCH_TIERS / GROMOS_BENCH_THREADS therefore narrow what actually
    # *runs*; the combinations left out are reported as skipped, and results
    # from a narrow run stay comparable with those from a full one.
    params = [list(TIERS), list(THREADS)]
    param_names = ["system", "threads"]

    # Generous: the large system on one thread, repeated, is slow.
    timeout = 7200.0

    def setup_cache(self):
        """Run every configuration once and return the parsed metrics.

        Returned value must be picklable -- it is a plain dict of floats keyed
        by ``(tier, threads)``. Configurations that do not apply to this build
        variant are simply absent, and the corresponding benchmarks then report
        as skipped.
        """
        var = _runner.variant()
        try:
            available = len(os.sched_getaffinity(0))
        except AttributeError:
            available = os.cpu_count() or 1

        # Fail loudly and early if the inputs are not the ones we expect: a
        # timeline built on drifting inputs measures nothing.
        _systems.load_manifest()

        results = {}
        for tier in _selected_tiers():
            for nthreads in _selected_threads():
                if not _supported(var, nthreads, available):
                    continue
                results["%s/%d" % (tier, nthreads)] = _runner.run_best(tier, nthreads,
                                                                      var=var)
        return results

    def setup(self, cache, system, threads):
        """Skip parameter combinations that do not apply to this build.

        NotImplementedError is asv's skip signal, but only when it comes from
        a setup function: asv catches it in do_setup, whereas the same
        exception raised inside a benchmark method is recorded as a failure.
        So the check has to live here rather than at the point of use, or the
        serial build reports a spurious failure for every thread count above
        one.
        """
        if "%s/%d" % (system, threads) not in cache:
            raise NotImplementedError(
                "%s at %d threads is not applicable to the %s build"
                % (system, threads, _runner.variant()))

    def _get(self, cache, system, threads, key):
        return cache["%s/%d" % (system, threads)].get(key)

    # -- headline -----------------------------------------------------------

    def track_sim_walltime(self, cache, system, threads):
        """Simulation wall time, excluding initialisation.

        Initialisation is excluded on purpose. For the large system, reading
        the 29 MB input coordinates takes roughly four times as long as the
        simulation being measured, so total wall time would be dominated by an
        I/O cost that has nothing to do with the integrator.
        """
        return self._get(cache, system, threads, "wall_simulation")

    track_sim_walltime.unit = "seconds"

    def track_ms_per_step(self, cache, system, threads):
        """Wall time per MD step -- comparable across systems of different size."""
        return self._get(cache, system, threads, "ms_per_step")

    track_ms_per_step.unit = "ms/step"

    # -- per-algorithm breakdown -------------------------------------------
    # These turn "the commit got slower" into "the pairlist got slower", which
    # is the difference between a graph you look at and one you can act on.

    def track_nonbonded(self, cache, system, threads):
        """Total time in the nonbonded interaction algorithm."""
        return self._get(cache, system, threads, "timer.NonBonded")

    track_nonbonded.unit = "seconds"

    def track_pairlist(self, cache, system, threads):
        """Pairlist construction, rebuilt every 5 steps in these setups."""
        return self._get(cache, system, threads, "subtimer.NonBonded.pairlist")

    track_pairlist.unit = "seconds"

    def track_nonbonded_shortrange(self, cache, system, threads):
        """Short-range solvent-solvent nonbonded work."""
        return self._get(cache, system, threads,
                         "subtimer.NonBonded.shortrange solvent-solvent")

    track_nonbonded_shortrange.unit = "seconds"

    def track_nonbonded_longrange(self, cache, system, threads):
        """Long-range solvent-solvent nonbonded work."""
        return self._get(cache, system, threads,
                         "subtimer.NonBonded.longrange solvent-solvent")

    track_nonbonded_longrange.unit = "seconds"

    def track_shake(self, cache, system, threads):
        """Constraint solving (SHAKE)."""
        return self._get(cache, system, threads, "timer.Shake")

    track_shake.unit = "seconds"

    # -- informational ------------------------------------------------------

    def track_ns_per_day(self, cache, system, threads):
        """Throughput. Higher is better, so ignore asv's regression flags here."""
        return self._get(cache, system, threads, "ns_per_day")

    track_ns_per_day.unit = "ns/day"
