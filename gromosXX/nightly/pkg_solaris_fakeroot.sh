source ${NIGHTHOME}/options.all
source ${NIGHTHOME}/options.${OPTIONNAME}

OS=`uname -s`"_"`uname -m`
source ${NIGHTHOME}/options.${OS}

SCRIPT=pkg_solaris
LOG=${NIGHT}/log/${SCRIPT}.log

LD_LIBRARY_PATH=/opt/csw/lib
export LD_LIBRARY_PATH

rm -rf ${TMP}/${NAME}_${SCRIPT}

ok=1

make DESTDIR=${TMP}/${NAME}_${SCRIPT} install >> ${LOG} 2>&1 || ok=0
cd ${TMP}/${NAME}_${SCRIPT} || ok=0
if [ ${ok} == 0 ] ; then
    echo "make install failed"
    echo "make install failed" >> ${LOG}
    exit 1
fi

(echo 'i pkginfo'; echo 'i copyright'; echo 'i depend'; pkgproto .) | sed 's: none : none /:' > prototype || ok=0
if [ ${ok} == 0 ] ; then
    echo "pkgproto failed"
    echo "pkgproto failed" >> ${LOG}
    exit 1
fi

unset LD_LIBRARY_PATH
