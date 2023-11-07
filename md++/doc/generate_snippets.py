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
# get lines containing the documentation example of a gromos block between 'INFO <<"' and '..\n";' 
# and output them as snippets that can be referred to in doxygen comments using @snippet

import sys, os

if len(sys.argv) < 2:
  sys.exit("USAGE: " + sys.argv[0] + " <filenames")

for f in sys.argv[1:]:
  with open(f, "r") as fi:
    snipping=0
    for line in fi:
      if snipping==0:
        if line.strip().replace(" ", "").replace("\t","").startswith('exampleblock<<"'):
          snipping=1
          blockname=line.split('"')[-2].split("\\n")[0]
          print("//! ["+blockname+"]")
          print(blockname)
      elif snipping==1:
        print(line.split('"')[-2].split("\\n")[0]+" ")
        if  line.strip().replace(" ", "").replace("\t","").startswith('exampleblock<<"END'):
          print("//! ["+blockname+"]\n")
          snipping=0
          
