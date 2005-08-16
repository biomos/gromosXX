#!/usr/bin/bash

mkdir -p /tmp/nightly_gromosXX_${BUILD}
cp ~/RELEASE/nightly/gromosXX-${VERSION}-${BUILD}.tar.gz /tmp/nightly_gromosXX_${BUILD}
cd /tmp/nightly_gromosXX_${BUILD}

${TAR} zxvf gromosXX-${VERSION}-${BUILD}.tar.gz > ${NIGHT}/log/distcheck_tar.log 2>&1
cd gromosXX-${VERSION}

./configure ${CONFIGURE_OPT} > ${NIGHT}/log/distcheck_configure.log 2>&1
make DISTCHECK_CONFIGURE_FLAGS="${CONFIGURE_OPT}" distcheck > ${NIGHT}/log/distcheck_make.log 2>&1

if [ $? == 0 ] ; then
    echo "distcheck succeeded" >> ${NIGHTLOG}
else
    echo "distcheck failed" >> ${NIGHTLOG}
fi

cd ~
rm -rf /tmp/nightly_gromosXX_${BUILD} > /dev/null 2>&1

