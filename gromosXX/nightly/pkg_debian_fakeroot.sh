
source ${NIGHTHOME}/options.all
OS=`uname -s`"_"`uname -m`
source ${NIGHTHOME}/options.${OS}

SCRIPT=pkg_solaris
LOG=${NIGHT}/log/${SCRIPT}.log

ok=1
rm -rf ${TMP}/${NAME}_pkg_debian
make DESTDIR=${TMP}/${NAME}_${SCRIPT} install >> ${LOG} 2>&1 || ok=0

cd ${TMP}/${NAME}_${SCRIPT} || ok=0
if [ ${ok} == 0 ] ; then
    echo "make install failed"
    echo "make install failed" >> ${LOG}
    exit 1
fi
