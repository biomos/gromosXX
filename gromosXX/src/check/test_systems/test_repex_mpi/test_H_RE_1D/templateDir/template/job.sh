gromosBIN=${3}
mpibin="" 


TOPO="input/omph_chcl3.top"
INPUTCRD="input/EQ_omph_chcl3_CD3OD-1_5.cnf"
IUNIT="input/HREMD_omph_chcl3_CD3OD-1_1.imd"
PTTOPO="input/dihedrals.pert"
OUTPUTCRD="test.cnf"
OUTPUTTRX="test_coord.trc"
OUTPUTTRE="test_ene.tre"
REPOUT="test_repout.out"
REPDAT="test_repdat.dat"
OUNIT="test.omd"

rm -f test*

echo -e "\t\t\tTesting: ${gromosBIN}"

${mpibin}mpirun -np ${2} ${gromosBIN}/repex_mpi  \
        @topo        ${TOPO} \
        @conf        ${INPUTCRD} \
        @input       ${IUNIT} \
        @pttopo      ${PTTOPO} \
        @fin         ${OUTPUTCRD} \
        @trc         ${OUTPUTTRX} \
        @tre         ${OUTPUTTRE}\
        @repout      ${REPOUT}\
        @repdat      ${REPDAT}\
	@verb re:8 \
        >            ${OUNIT}


