#!/usr/bin/env python3
"""Convert relative #include "..." to rooted #include <...> anchored at src/."""
import os
import re
import sys

MD_ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))
SRC_ROOT = os.path.join(MD_ROOT, "src")

# Mirrors the -I order used by both build systems
INCLUDE_ROOTS = [SRC_ROOT, MD_ROOT]

EXTENSIONS = {".cc", ".cpp", ".h", ".cu", ".cuh"}
INCLUDE_RE = re.compile(r'^(\s*#\s*include\s*)"([^"]+)"')


def resolve_include(file_dir, rel):
    """Return path relative to SRC_ROOT, or None if not resolvable."""
    # 1. Try relative to the including file's directory
    candidate = os.path.realpath(os.path.join(file_dir, rel))
    if candidate.startswith(SRC_ROOT + os.sep) and os.path.exists(candidate):
        return os.path.relpath(candidate, SRC_ROOT)

    # 2. Resolved path is outside src/ (e.g. excess ../ in qmmm files).
    #    Fall back to searching the include roots for the trailing path
    #    component (everything after the last ../ sequence).
    parts = rel.replace("\\", "/").split("/")
    i = 0
    while i < len(parts) and parts[i] == "..":
        i += 1
    trailing = "/".join(parts[i:])
    for root in INCLUDE_ROOTS:
        full = os.path.join(root, trailing)
        if os.path.exists(full):
            return os.path.relpath(os.path.realpath(full), SRC_ROOT)

    return None


def convert_file(filepath):
    file_dir = os.path.dirname(os.path.realpath(filepath))
    with open(filepath) as f:
        lines = f.readlines()
    changed = False
    result = []
    for line in lines:
        m = INCLUDE_RE.match(line)
        if m:
            prefix, rel = m.group(1), m.group(2)
            rooted = resolve_include(file_dir, rel)
            if rooted is not None:
                new_line = f'{prefix}<{rooted}>\n'
                if new_line != line:
                    line = new_line
                    changed = True
        result.append(line)
    if changed:
        with open(filepath, "w") as f:
            f.writelines(result)
        print(f"Updated: {os.path.relpath(filepath, SRC_ROOT)}")


roots = [SRC_ROOT]
if len(sys.argv) > 1:
    roots = [os.path.realpath(a) for a in sys.argv[1:]]

for root in roots:
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if os.path.splitext(fname)[1] in EXTENSIONS:
                convert_file(os.path.join(dirpath, fname))
