"""Copy the benchmark systems into a local cache.

This is the only part of the suite that reads the upstream protein set. The
upstream location is shared storage that may be unavailable, slow, or mutable;
benchmark inputs must be none of those. After staging, every other module reads
exclusively from the local cache.

The manifest records a sha256 per file. The harness verifies it before every
run, because benchmark inputs must be byte-identical across every commit in the
timeline -- if the inputs drift, the timeline stops measuring the code.
"""

import argparse
import hashlib
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from benchmarks import _systems


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return "%.1f %s" % (n, unit)
        n /= 1024.0


def stage(tiers, force):
    dest_root = _systems.cache_root()
    manifest = {"source": _systems.source_root(), "systems": {}}

    manifest_path = os.path.join(dest_root, _systems.MANIFEST_NAME)
    if os.path.isfile(manifest_path) and not force:
        with open(manifest_path) as fh:
            manifest = json.load(fh)
        manifest.setdefault("systems", {})

    total = 0
    for tier in tiers:
        spec = _systems.SYSTEMS[tier]
        pdb = spec["pdb"]
        src = _systems.source_files(pdb)
        dst = _systems.staged_files(pdb)

        missing = [p for p in src.values() if not os.path.isfile(p)]
        if missing:
            raise SystemExit(
                "cannot stage %r (%s) -- missing upstream files:\n  %s"
                % (tier, pdb, "\n  ".join(missing)))

        os.makedirs(os.path.dirname(dst["top"]), exist_ok=True)
        entry = {"pdb": pdb, "tier": tier, "atoms": spec["atoms"], "files": {}}

        for kind in ("top", "cnf", "imd"):
            s, d = src[kind], dst[kind]
            if force or not os.path.isfile(d) or os.path.getsize(d) != os.path.getsize(s):
                print("  copying %-4s %s" % (kind, os.path.basename(s)), flush=True)
                shutil.copyfile(s, d)
            size = os.path.getsize(d)
            total += size
            entry["files"][kind] = {
                "name": os.path.basename(d),
                "size": size,
                "sha256": sha256(d),
            }

        manifest["systems"][tier] = entry
        print("staged %-6s %-5s %7d atoms" % (tier, pdb, spec["atoms"]), flush=True)

    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print("\ncache:    %s" % dest_root)
    print("manifest: %s" % manifest_path)
    print("total:    %s" % human(total))


def verify(tiers):
    """Re-hash the staged files and compare against the manifest."""
    manifest = _systems.load_manifest()
    bad = []
    for tier in tiers:
        entry = manifest["systems"].get(tier)
        if entry is None:
            bad.append("%s: not in manifest" % tier)
            continue
        files = _systems.staged_files(entry["pdb"])
        for kind, rec in entry["files"].items():
            path = files[kind]
            if not os.path.isfile(path):
                bad.append("%s/%s: missing" % (tier, kind))
            elif sha256(path) != rec["sha256"]:
                bad.append("%s/%s: checksum mismatch" % (tier, kind))
        print("verified %-6s %s" % (tier, entry["pdb"]), flush=True)
    if bad:
        raise SystemExit("staged systems do not match the manifest:\n  "
                         + "\n  ".join(bad))
    print("\nall staged systems match the manifest")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--systems", default="all",
                    help="comma-separated tiers to stage (default: all)")
    ap.add_argument("--force", action="store_true",
                    help="re-copy even if the destination already looks correct")
    ap.add_argument("--verify", action="store_true",
                    help="only re-check staged files against the manifest")
    args = ap.parse_args()

    if args.systems == "all":
        tiers = list(_systems.TIERS)
    else:
        tiers = [t.strip() for t in args.systems.split(",") if t.strip()]
        unknown = [t for t in tiers if t not in _systems.SYSTEMS]
        if unknown:
            raise SystemExit("unknown tier(s): %s" % ", ".join(unknown))

    if args.verify:
        verify(tiers)
    else:
        stage(tiers, args.force)


if __name__ == "__main__":
    main()
