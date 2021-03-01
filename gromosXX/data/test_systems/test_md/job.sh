
rm test*;

gromosDIR=/cluster/work/igc/bschroed/gromos/gromos_compiled/gromosXX_more_Parallel_newEnvMPI2/bin
${gromosDIR}/md \
    @topo ./input/peptide_2Cl_54a7.top\
    @conf ./input/eq_peptide_5.cnf\
    @input ./input/md.imd \
    @fin test.cnf \
    @trc test.trc \
    @tre test.tre \
    @trs test.trs \
    @verb 0 \
    ${1} >test.omd

#    @verb 5 \
