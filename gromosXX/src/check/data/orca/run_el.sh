#!/bin/bash

/home/fpultar/src/gromos-standard/gromosXX/gromosXX/build-vscode/program/md \
@topo md.top \
@conf md.cnf \
@input md_el.imd \
@qmmm md_el.qmmm \
@fin md_final.cnf \
@trc md_el.trc \
@tre md_el.tre > md_el.omd
