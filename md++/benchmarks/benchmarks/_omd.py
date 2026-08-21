"""Parse timing metrics out of an md++ output (omd) file.

md++ already reports everything a performance benchmark needs, which is why the
suite runs each simulation once and scrapes many metrics rather than timing the
process externally. Two sources, both written by ``program/md.cc``:

* a summary block::

      Wall time initialisation (s):        0.660
      Wall time simulation (s):            1.619
      Wall time total (s):                 2.279
      Performance (ns/day):               21.341

* a ``TIMING`` block from ``util::Algorithm_Timer::print``
  (``src/util/timing.cc``) giving a per-algorithm breakdown, where algorithms
  with sub-timers list them indented underneath::

      TIMING
                NonBonded                              0.031     100.00%
                   - pairlist                                    0.009  28.22%
                   - shortrange solvent-solvent                  0.008  27.58%
      END

Sub-timer names contain spaces ("longrange solvent-solvent"), so the name is
matched non-greedily up to the numeric columns rather than by splitting.
"""

import re

_WALL = re.compile(r"^Wall time (\w+) \(s\):\s+([0-9.eE+-]+)\s*$")
_PERF = re.compile(r"^Performance \(ns/day\):\s+([0-9.eE+-]+)\s*$")

# "             - <name padded to 35> <value> <percent>%"
_SUBTIMER = re.compile(r"^\s+-\s+(\S.*?)\s+([0-9.]+)\s+([0-9.]+)%\s*(\[.*\])?\s*$")

# "          <name padded to 30> <value> [<percent>%]"
_MAINTIMER = re.compile(r"^\s{6,}([A-Za-z]\S*)\s+([0-9.]+)\s*(?:([0-9.]+)%)?\s*$")

SUCCESS_MARKER = "MD++ finished successfully"


class OmdError(RuntimeError):
    """The run did not finish successfully, or produced no usable timings."""


def parse(text):
    """Return a metrics dict parsed from omd text.

    Keys:
      ``wall_initialisation`` / ``wall_simulation`` / ``wall_total`` -- seconds
      ``ns_per_day``                                                -- ns/day
      ``timer.<Algorithm>``                                         -- seconds
      ``subtimer.<Algorithm>.<sub name>``                           -- seconds
    """
    if SUCCESS_MARKER not in text:
        raise OmdError(_failure_detail(text))

    metrics = {}
    in_timing = False
    current = None

    for line in text.splitlines():
        if line.strip() == "TIMING":
            in_timing = True
            continue
        if in_timing:
            if line.strip() == "END":
                in_timing = False
                current = None
                continue
            m = _SUBTIMER.match(line)
            if m and current:
                name = m.group(1).strip()
                # "unaccounted time" is bookkeeping, not an algorithm.
                if name != "unaccounted time":
                    metrics["subtimer.%s.%s" % (current, name)] = float(m.group(2))
                continue
            m = _MAINTIMER.match(line)
            if m:
                current = m.group(1)
                metrics["timer.%s" % current] = float(m.group(2))
            continue

        m = _WALL.match(line)
        if m:
            metrics["wall_%s" % m.group(1)] = float(m.group(2))
            continue
        m = _PERF.match(line)
        if m:
            metrics["ns_per_day"] = float(m.group(1))

    if "wall_simulation" not in metrics:
        raise OmdError("no 'Wall time simulation' line in omd output; "
                       "the run reported success but produced no timings")
    return metrics


def parse_file(path):
    with open(path, errors="replace") as fh:
        return parse(fh.read())


def _failure_detail(text):
    """Build a useful message from md++'s own error reporting."""
    errors = [l.strip() for l in text.splitlines()
              if re.match(r"^\s*(ERROR|Errors during)", l)]
    detail = "\n  ".join(errors[:8]) if errors else "(no ERROR lines found)"
    return ("md++ did not report %r -- the run failed.\n  %s"
            % (SUCCESS_MARKER, detail))
