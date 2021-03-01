out_prefix="${1}"
gromosDIR=/cluster/work/igc/bschroed/gromos/gromos_compiled/gromosXX_more_Parallel_newEnvMPI2/bin

# User specific aliases and functions
. ~/new_env.sh



mpirun ${gormosDIR}/md_mpi \
    @topo ./input/peptide_2Cl_54a7.top\
    @conf ./input/eq_peptide_5.cnf\
    @input ./input/md.imd \
    @fin test.cnf \
    @trc test.trc \
    @tre test.tre \
    @trs test.trs \
     >test.omd

#    @verb 5 \
