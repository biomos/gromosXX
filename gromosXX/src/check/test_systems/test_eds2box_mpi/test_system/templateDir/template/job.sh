 #!/bin/bash

. ../../../gromos_settings.sh

# first we set some variables
NAME="TheGreatTester"

#remove previous files
rm -f test*;

SIMULDIR=${PWD}
gromDIR="${3}"


# set the input files for box 1
TOPO1=${SIMULDIR}/input/topo/1000_SPC.top
IUNIT1=${SIMULDIR}/input/imd_eds_1.inp1
INPUTCRD1=${SIMULDIR}/input/md_eds_0.fin1
PTTOPO1=${SIMULDIR}/input/topo/pert_eds.ptp1
# set the input files for box 2
TOPO2=${SIMULDIR}/input/topo/900_SPC_100_MeOH.top
IUNIT2=${SIMULDIR}/input/imd_eds_1.inp2
INPUTCRD2=${SIMULDIR}/input/md_eds_0.fin2
PTTOPO2=${SIMULDIR}/input/topo/pert_eds.ptp2


#set the output files for box 1
outprefix=${1}
OUNIT1=${outprefix}_1.out
OUTPUTCRD1=${outprefix}_1.fin1
OUTPUTTRX1=${outprefix}_1.trj1
OUTPUTTRE1=${outprefix}_1.tre1
#set the output files for box 2
#OUNIT2=md_eds_1.out2
OUTPUTCRD2=${outprefix}_1.fin2
OUTPUTTRX2=${outprefix}_1.trj2
OUTPUTTRE2=${outprefix}_1.tre2


if [ ${gromDIR} == *"debug"* ]
then
mpirun -np ${2}  ${gromDIR}/eds_2box_mpi \
        @topo1        ${TOPO1} \
        @topo2        ${TOPO2} \
        @conf1        ${INPUTCRD1} \
        @conf2        ${INPUTCRD2} \
        @input1       ${IUNIT1} \
        @input2       ${IUNIT2} \
        @pttopo1      ${PTTOPO1} \
        @pttopo2      ${PTTOPO2} \
        @fin1         ${OUTPUTCRD1} \
        @fin2         ${OUTPUTCRD2} \
        @trc1         ${OUTPUTTRX1} \
        @trc2         ${OUTPUTTRX2} \
        @tre1         ${OUTPUTTRE1}\
        @tre2         ${OUTPUTTRE2}\
        @verb   2
        >             ${OUNIT1}   || exit 1
else

mpirun -np  ${2}  ${gromDIR}/eds_2box_mpi \
        @topo1        ${TOPO1} \
        @topo2        ${TOPO2} \
        @conf1        ${INPUTCRD1} \
        @conf2        ${INPUTCRD2} \
        @input1       ${IUNIT1} \
        @input2       ${IUNIT2} \
        @pttopo1      ${PTTOPO1} \
        @pttopo2      ${PTTOPO2} \
        @fin1         ${OUTPUTCRD1} \
        @fin2         ${OUTPUTCRD2} \
        @trc1         ${OUTPUTTRX1} \
        @trc2         ${OUTPUTTRX2} \
        @tre1         ${OUTPUTTRE1}\
        @tre2         ${OUTPUTTRE2}\
        >             ${OUNIT1}   || exit 1
fi


echo done

