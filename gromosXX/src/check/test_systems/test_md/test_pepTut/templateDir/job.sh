#Load Gromos Settings
. ../../../gromos_settings.sh

#remove previous files
rm -f test*;
gromDIR=${1} 

if [ ${gromDIR} == *debug* ]
then
    ${gromDIR}/md \
        @topo ./input/peptide_2Cl_54a7.top\
        @conf ./input/eq_peptide_5.cnf\
        @input ./input/md.imd \
        @fin test.cnf \
        @trc test.trc \
        @tre test.tre \
        @trs test.trs \
        @verb 5 \
        ${1} >test.omd || exit 1
else

    ${gromDIR}/md \
        @topo ./input/peptide_2Cl_54a7.top\
        @conf ./input/eq_peptide_5.cnf\
        @input ./input/md.imd \
        @fin test.cnf \
        @trc test.trc \
        @tre test.tre \
        @trs test.trs \
        ${1} >test.omd || exit 1
fi
