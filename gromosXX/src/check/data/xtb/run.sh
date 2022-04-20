#!/bin/bash
source ~/opt/gcc-8.2.0/activate

/home/fpultar/src/gromos-standard/gromosXX/gromosXX/build-vscode/program/md \
@topo md.top \
@conf md.cnf \
@input md.imd \
@qmmm md.qmmm \
@fin md_final.cnf \
@trc mc.trc \
@tre md.tre > md.omd
