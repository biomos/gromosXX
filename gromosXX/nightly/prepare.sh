#/bin/sh

echo "build ${BUILD}"
echo "build ${BUILD}" >> ${NIGHTLOG}

echo `date "+%d.%m.%y %T"`"     distribution creation started" >> ${NIGHTLOG}

BUILDDIR=${TMP}/nightly_${NAME}_${BUILD}
SCRIPT=prepare
LOG=${NIGHT}/log/${SCRIPT}.log

cat /dev/null > ${LOG}

rm -rf ${BUILDDIR}
mkdir -p ${BUILDDIR}

cd ${BUILDDIR}

ok=1
echo "     cvs co ${NAME}" >> ${LOG}
cvs co ${NAME} >> /dev/null 2>&1 || ok=0

if [ ${ok} == 0 ] ; then
    echo "cvs co ${NAME} failed"
    echo "cvs co ${NAME} failed" >> ${LOG}
    exit 1
fi

cd ${NAME}
if [ -x Config.sh ] ; then
    ./Config.sh >> ${LOG} 2>&1 || ok=0

    if [ ${ok} == 0 ] ; then
	echo "Config.sh failed!"
	echo "Config.sh failed!" >> ${LOG}
	exit 1
    fi
fi

./configure ${CONFIGURE_OPT} >> ${LOG} 2>&1 || ok=0
if [ ${ok} == 0 ] ; then
    echo "./configure ${CONFIGURE_OPT} failed"
    echo "./configure ${CONFIGURE_OPT} failed" >> ${LOG}
    exit 1
fi

make dist >> ${LOG} 2>&1 || ok=0
if [ ${ok} == 0 ] ; then
    echo "make dist failed"
    echo "make dist failed" >> ${LOG}
    exit 1
fi

echo `date "+%d.%m.%y %T"`"     distribution created" >> ${NIGHTLOG}

cp ${NAME}-${VERSION}.tar.gz ${NIGHT}/${NAME}-${VERSION}-${BUILD}.tar.gz

cd ~
rm -rf ${BUILDDIR}
