#!/bin/sh

# first we set some variables
NAME=`TheGreatTester`

PROGRAM=eds_2box
SIMULDIR=${PWD}

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
OUNIT1=tmd_eds_1.out
OUTPUTCRD1=tmd_eds_1.fin1
OUTPUTTRX1=tmd_eds_1.trj1
OUTPUTTRE1=tmd_eds_1.tre1
#set the output files for box 2
#OUNIT2=md_eds_1.out2
OUTPUTCRD2=tmd_eds_1.fin2
OUTPUTTRX2=tmd_eds_1.trj2
OUTPUTTRE2=tmd_eds_1.tre2


MDOK=1

echo "START THE FUN!"
#${PROGRAM} @verb algorithm:integration:7 \
#mpirun ${PROGRAM} \
mpirun -np 2 ${PROGRAM} \
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
        >             ${OUNIT1}
grep "finished successfully" ${OUNIT1} > /dev/null || MDOK=0

uname -a >> ${OUNIT1}

# compress some files
gzip ${OUTPUTTRX1}
gzip ${OUTPUTTRE1}
gzip ${OUTPUTTRX2}
gzip ${OUTPUTTRE2}


# clean up after us
if `test ${OK} -eq 0`; then
  uname -a > mess;
  echo 'cp failed for WT_mp, run 1' >> mess;
  Mail -s "ERROR" ${NAME} < mess;
  cd ${SIMULDIR};
else
  cd ${SIMULDIR};
  rm ${WORKDIR}/*;
  rmdir ${WORKDIR};
fi

echo done

