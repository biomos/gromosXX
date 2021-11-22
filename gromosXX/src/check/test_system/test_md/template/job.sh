

rm test*;
gromDIR=${1} 

${gromDIR}/md \
    @topo ./input/peptide_2Cl_54a7.top\
    @conf ./input/eq_peptide_5.cnf\
    @input ./input/md.imd \
    @fin test.cnf \
    @trc test.trc \
    @tre test.tre \
    @trs test.trs \
    @verb 5 \
    ${1} >test.omd

#    @verb 5 \
