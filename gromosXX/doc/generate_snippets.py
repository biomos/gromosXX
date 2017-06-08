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
          
