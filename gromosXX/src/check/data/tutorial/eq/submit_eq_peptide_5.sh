#!/bin/bash
export SIMULDIR=$(pwd);
bsub -n 16 -W 24:00 -R "rusage[scratch=1000]" -J eq_peptide_5 -e ${SIMULDIR}/eq_peptide_5.err -o ${SIMULDIR}/eq_peptide_5.out < eq_peptide_5.run
