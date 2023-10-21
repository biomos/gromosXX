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
''' convert old gromos doxygen block descriptions 
    to new ones (INFO<<"";)
'''
import sys, os, glob

if len(sys.argv)<2:
  sys.exit("USAGE: "+sys.argv[0] + " <source file>")

intro='''
  std::stringstream INFO;
  // lines starting with 'INFO<<"' and ending with '\\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag\n'''

with open(sys.argv[1]) as f: 
  read=False
  for line in f:
    if line.find("@verbatim") != -1:
      read=True
      infotxt=(intro)
      line=next(f)
      blockname=line.strip()
      print(" * @snippet snippets/snippets.cc "+blockname)
      infotxt+='  INFO << "' +line.rstrip()+'\\n";\n'
    elif line.find("@endverbatim")!= -1:
      read=False
   # elif line.find("DEBUG(8,") != -1 and line.find("reading") != -1:
    elif line.find("DEBUG(10,") != -1 and line.find("block") != -1:
      print(line.rstrip())
      print(infotxt)
    elif read and not line.strip()=="":
      infotxt+='  INFO << "' +line.rstrip()+'\\n";\n'      
    else:
      print(line.rstrip())

   
