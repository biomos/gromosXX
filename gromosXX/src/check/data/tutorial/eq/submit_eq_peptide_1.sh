#!/bin/bash
export SIMULDIR=$(pwd);
bsub -n 16 -W 24:00 -R "rusage[scratch=1000]" -J eq_peptide_1 -e ${SIMULDIR}/eq_peptide_1.err -o ${SIMULDIR}/eq_peptide_1.out < eq_peptide_1.run
