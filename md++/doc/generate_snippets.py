"""
This file is part of GROMOS.

Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
See <https://www.gromos.net> for details.

GROMOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

#! /usr/bin/env python
# Extract lines containing the documentation example of a GROMOS block
# between 'exampleblock <<"' and '..\n";' and output them as snippets
# that can be referred to in Doxygen comments using @snippet.

import sys
import os

if len(sys.argv) < 3:
    sys.exit("USAGE: " + sys.argv[0] + " <output_dir> <filenames>")

output_dir = sys.argv[1]
input_files = sys.argv[2:]


print(output_dir)
print(input_files)
# Ensure the output directory exists and is valid
if os.path.isfile(output_dir):
    sys.exit(f"ERROR: Output directory '{output_dir}' conflicts with an existing file.")
os.makedirs(output_dir, exist_ok=True)

for f in input_files:
    if not os.path.isfile(f):
        print(f"WARNING: Input file '{f}' does not exist. Skipping.")
        continue

    with open(f, "r") as fi:
        snipping = 0
        blockname = None
        output_lines = []
        for line in fi:
            if snipping == 0:
                if line.strip().replace(" ", "").replace("\t", "").startswith('exampleblock<<"'):
                    snipping = 1
                    blockname = line.split('"')[-2].split("\\n")[0]
                    output_lines.append(f"//! [{blockname}]")
            elif snipping == 1:
                output_lines.append(line.split('"')[-2].split("\\n")[0] + " ")
                if line.strip().replace(" ", "").replace("\t", "").startswith('exampleblock<<"END'):
                    output_lines.append(f"//! [{blockname}]\n")
                    snipping = 0
                    # Write the snippet to a file
                    snippet_file = os.path.join(output_dir, f"{blockname}.snippet")
                    with open(snippet_file, "w") as fo:
                        fo.write("\n".join(output_lines))
                    output_lines = []