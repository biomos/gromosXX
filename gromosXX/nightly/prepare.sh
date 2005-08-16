#/bin/sh

echo "build ${BUILD}"

mkdir -p /tmp/nightly_gromosXX_${BUILD}
cd /tmp/nightly_gromosXX_${BUILD}

CVSROOT=:pserver:igc@igc.ethz.ch:/home/cvs/gromosXX
cvs co gromosXX

cd gromosXX
./Config.sh
./configure ${CONFIGURE_OPT}

make dist && echo "gromosXX-${VERSION}.tar.gz created" >> ${NIGHTLOG}

cp gromosXX-${VERSION}.tar.gz ~/RELEASE/nightly/gromosXX-${VERSION}-${BUILD}.tar.gz
cd ~
rm -rf /tmp/nightly_gromosXX_${BUILD}
