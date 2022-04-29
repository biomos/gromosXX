#!/bin/bash

/home/fpultar/src/gromos-standard/gromosXX/gromosXX/build-vscode/program/md \
@topo md.top \
@conf md.cnf \
@input md_mc.imd \
@qmmm md_mc.qmmm \
@fin md_final.cnf \
@trc md_mc.trc \
@tre md_mc.tre > md_mc.omd
