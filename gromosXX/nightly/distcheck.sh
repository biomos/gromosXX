#!/usr/bin/bash

SCRIPT=distcheck
LOG=${NIGHT}/log/${SCRIPT}.log
BUILDDIR=${TMP}/nightly_${NAME}_${BUILD}
INSTALLDIR=${TMP}/${NAME}_${SCRIPT}

cat /dev/null > ${LOG}

echo `date "+%d.%m.%y %T"`"     distcheck started" >> ${NIGHTLOG}

rm -rf ${BUILDDIR}
mkdir -p ${BUILDDIR}

ok=1
cp ${NIGHT}/${NAME}-${VERSION}-${BUILD}.tar.gz ${BUILDDIR} || ok=0
cd ${BUILDDIR} || ok=0

if [ ${ok} == 0 ] ; then
    echo "preparing directory failed"
    echo "preparing directory failed" >> ${LOG}
    exit 1
fi

${TAR} zxvf ${NAME}-${VERSION}-${BUILD}.tar.gz >> ${LOG} 2>&1 || ok=0
cd ${NAME}-${VERSION} || ok=0
if [ ${ok} == 0 ] ; then
    echo "unpacking sources failed"
    echo "unpacking sources failed" >> ${LOG}
    exit 1
fi

./configure ${CONFIGURE_OPT} >> ${LOG} 2>&1 || ok=0
if [ ${ok} == 0 ] ; then
    echo "configure ${CONFIGURE_OPT} failed"
    echo "configure ${CONFIGURE_OPT} failed" >> ${LOG}
    exit 1
fi

make DISTCHECK_CONFIGURE_FLAGS="${CONFIGURE_OPT}" distcheck >> ${LOG} 2>&1 || ok=0
if [ ${ok} == 0 ] ; then
    echo "make DISTCHECK_CONFIGURE_FLAGS=${CONFIGURE_OPT} distcheck failed"
    echo "make DISTCHECK_CONFIGURE_FLAGS=${CONFIGURE_OPT} distcheck failed" >> ${LOG}
    exit 1
fi

echo `date "+%d.%m.%y %T"`"     distcheck succeeded" >> ${NIGHTLOG}

cd ~
rm -rf ${BUILDDIR}

