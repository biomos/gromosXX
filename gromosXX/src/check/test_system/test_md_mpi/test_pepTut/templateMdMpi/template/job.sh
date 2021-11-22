out_prefix="${1}"
gromosBIN="${3}"



mpirun -n ${2} ${gromosBIN}/md_mpi \
    @topo ./input/peptide_2Cl_54a7.top\
    @conf ./input/eq_peptide_5.cnf\
    @input ./input/md.imd \
    @fin test.cnf \
    @trc test.trc \
    @tre test.tre \
    @trs test.trs \
     >test.omd || exit 1

#    @verb 5 \
