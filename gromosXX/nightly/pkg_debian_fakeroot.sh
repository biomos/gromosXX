
source ${NIGHTHOME}/options.all
source ${NIGHTHOME}/options.${OPTIONNAME}
OS=`uname -s`"_"`uname -m`
source ${NIGHTHOME}/options.${OS}

SCRIPT=pkg_debian
LOG=${NIGHT}/log/${SCRIPT}.log
INSTALLDIR=${TMP}/${NAME}_${SCRIPT}

ok=1
rm -rf ${INSTALLDIR}
make DESTDIR=${INSTALLDIR} install >> ${LOG} 2>&1 || ok=0

cd ${INSTALLDIR} || ok=0
if [ ${ok} == 0 ] ; then
    echo "make install failed"
    echo "make install failed" >> ${LOG}
    exit 1
fi
