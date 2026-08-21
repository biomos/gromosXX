"""Minimal block-aware reader/writer for GROMOS imd parameter files.

An imd file is a flat sequence of ``NAME ... END`` blocks. Within a block,
lines starting with ``#`` are comments describing the parameter columns; the
rest are whitespace-separated values.

Editing has to be block-aware rather than line-based, because some blocks take
a *variable* number of parameters. ``INNERLOOP`` is the one that matters here:
with ``NTILM=0`` only two values may appear, and supplying the four that a CUDA
run needs makes md++ fail during initialisation with

    ERROR In_Parameter : Left-over parameters in INNERLOOP block

so a naive regex substitution that preserved the original column count would
silently produce an unrunnable input for CPU runs.
"""

import re

_BLOCK_START = re.compile(r"^([A-Z][A-Z0-9_]*)\s*$")


class ImdFile(object):
    """An imd file as an ordered list of blocks, preserving comments."""

    def __init__(self, text):
        self.blocks = []          # list of (name, list-of-lines) in file order
        self._parse(text)

    @classmethod
    def load(cls, path):
        with open(path) as fh:
            return cls(fh.read())

    def _parse(self, text):
        name = None
        body = []
        preamble = []
        for line in text.splitlines():
            if name is None:
                m = _BLOCK_START.match(line)
                if m:
                    name = m.group(1)
                    body = []
                else:
                    # Blank lines and stray text between blocks.
                    preamble.append(line)
                continue
            if line.strip() == "END":
                self.blocks.append((name, body))
                name = None
            else:
                body.append(line)
        if name is not None:
            raise ValueError("unterminated block %r in imd file" % name)
        self._trailing = preamble

    def has(self, name):
        return any(n == name for n, _ in self.blocks)

    def comments(self, name):
        """The ``#`` lines of a block, which document its parameter columns."""
        for n, body in self.blocks:
            if n == name:
                return [l for l in body if l.lstrip().startswith("#")]
        raise KeyError(name)

    def values(self, name):
        """Whitespace-separated non-comment tokens of a block."""
        for n, body in self.blocks:
            if n == name:
                toks = []
                for line in body:
                    if not line.lstrip().startswith("#"):
                        toks.extend(line.split())
                return toks
        raise KeyError(name)

    def set_values(self, name, values):
        """Replace a block's parameters, keeping its comment lines intact.

        ``values`` fully determines the new parameter count, which is what
        makes variable-arity blocks work.
        """
        payload = "    " + "  ".join(str(v) for v in values)
        for i, (n, body) in enumerate(self.blocks):
            if n == name:
                kept = [l for l in body if l.lstrip().startswith("#")]
                self.blocks[i] = (n, kept + [payload])
                return
        raise KeyError("no %s block in imd file" % name)

    def dumps(self):
        out = []
        for name, body in self.blocks:
            out.append(name)
            out.extend(body)
            out.append("END")
        return "\n".join(out) + "\n"

    def save(self, path):
        with open(path, "w") as fh:
            fh.write(self.dumps())


def prepare_benchmark_imd(src, dest, nstlim, innerloop):
    """Derive a benchmark imd from a production one.

    Three edits, each load-bearing:

    * ``STEP``      -- set NSTLIM so the run lands in the target duration, and
                       force the start time to 0 so runs are comparable.
    * ``WRITETRAJ`` -- zero every output frequency. Left alone, a production imd
                       writes coordinate and energy trajectories and the
                       benchmark measures the filesystem as much as the physics.
    * ``INNERLOOP`` -- select the CPU or CUDA nonbonded kernel. This is a
                       *runtime* choice, so one CUDA-capable build serves both.

    Returns the resulting :class:`ImdFile`.
    """
    imd = ImdFile.load(src)

    step = imd.values("STEP")
    if len(step) != 3:
        raise ValueError("unexpected STEP block arity %d (expected NSTLIM T DT)"
                         % len(step))
    imd.set_values("STEP", [nstlim, 0, step[2]])          # keep the timestep

    if imd.has("WRITETRAJ"):
        imd.set_values("WRITETRAJ", [0] * len(imd.values("WRITETRAJ")))

    imd.set_values("INNERLOOP", innerloop)

    imd.save(dest)
    return imd
