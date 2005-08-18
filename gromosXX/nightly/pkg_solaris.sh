#!/usr/bin/bash

SCRIPT=pkg_solaris
LOG=${NIGHT}/log/${NAME}_${SCRIPT}.log
BUILDDIR=${TMP}/nightly_${NAME}_${BUILD}
INSTALLDIR=${TMP}/${NAME}_${SCRIPT}

cat /dev/null > ${LOG}

echo `date "+%d.%m.%y %T"`"   ${NAME} solaris package creation started" >> ${NIGHTLOG}

echo "creating a solaris package" >> ${LOG}
echo "     preparing directory" >> ${LOG}
ok=1
rm -rf ${BUILDDIR}
mkdir -p ${BUILDDIR} || ok=0
cp ${NIGHT}/${NAME}-${VERSION}-${BUILD}.tar.gz ${BUILDDIR} || ok=0
cd ${BUILDDIR} || ok=0

if [ ${ok} == 0 ] ; then
    echo "preparing directory failed"
    echo "preparing directory failed" >> ${LOG}
    exit 1
fi

echo "    tar" >> ${LOG}
${TAR} zxvf ${NAME}-${VERSION}-${BUILD}.tar.gz > /dev/null 2>&1 || ok=0
cd ${NAME}-${VERSION} || ok=0

if [ ${ok} == 0 ] ; then
    echo "${TAR} zxvf ${NAME}-${VERSION}-${BUILD}.tar.gz failed"
    echo "${TAR} zxvf ${NAME}-${VERSION}-${BUILD}.tar.gz failed" >> ${LOG}
    exit 1
fi

echo "" >> ${LOG}
echo "     configuring" >> ${LOG}
./configure ${CONFIGURE_REL_OPT} >> ${LOG} 2>&1 || ok=0
echo "     make" >> ${LOG}
make >> ${LOG} 2>&1 || ok=0

if [ ${ok} == 0 ] ; then
    echo "configure or make failed"
    echo "configure or make failed" >> ${LOG}
    cp config.log ${NIGHT}/log/${SCRIPT}.config.log
    exit 1
fi

echo "     enter fakeroot" >> ${LOG}
fakeroot ${BASH} ${NIGHTHOME}/${SCRIPT}_fakeroot.sh || ok=0

if [ ${ok} == 0 ] ; then
    echo "enter fakeroot failed"
    echo "enter fakeroot failed" >> ${LOG}
    exit 1
fi

cd ${INSTALLDIR} || ok=0
if [ ${ok} == 0 ] ; then
    echo "cd ${INSTALLDIR} failed"
    echo "cd ${INSTALLDIR} failed" >> ${LOG}
    exit 1
fi

cp ${BUILDDIR}/${NAME}-${VERSION}/packaging/solaris.depend depend
cp ${BUILDDIR}/${NAME}-${VERSION}/packaging/copyright .

ARCH=`uname -p`

echo "PKG=\"${PKGNAME}\"" > pkginfo || ok=0
echo "NAME=\"${DESCRIPTION}\"" >> pkginfo || ok=0
echo "VERSION=\"${VERSION} ${BUILD}\"" >> pkginfo || ok=0
echo "ARCH=\"${ARCH}\"" >> pkginfo || ok=0
echo "CLASSES=\"none\"" >> pkginfo || ok=0
echo "CATEGORY=\"develop\"" >> pkginfo || ok=0
echo "VENDOR=\"${MAINTAINER}\"" >> pkginfo || ok=0
echo "PSTAMP=\"${DATE}\"" >> pkginfo || ok=0
echo "BASEDIR=/" >> pkginfo || ok=0

if [ ${ok} == 0 ] ; then
    echo "creating pkginfo failed"
    echo "creating pkginfo failed" >> ${LOG}
    exit 1
fi

pkgmk -r ${TMP}/${NAME}_${SCRIPT} >> ${LOG} 2>&1 || ok=0
if [ ${ok} == 0 ] ; then
    echo "pkgmk -r ${TMP}/${NAME}_${SCRIPT} failed"
    echo "pkgmk -r ${TMP}/${NAME}_${SCRIPT} failed" >> ${LOG}
    exit 1
fi

PKGLONGNAME=${PKGNAME}-${VERSION}-${BUILD}-${OS}.pkg
pkgtrans -s /var/spool/pkg ${TMP}/${PKGLONGNAME} ${PKGNAME} >> ${LOG} || ok=0
if [ ${ok} == 0 ] ; then
    echo "pkgtrans -s /var/spool/pkg ${TMP}/${PKGLONGNAME}.pkg ${PKGNAME} failed"
    echo "pkgtrans -s /var/spool/pkg ${TMP}/${PKGLONGNAME}.pkg ${PKGNAME} failed" >> ${LOG}
    exit 1
fi

mv ${TMP}/${PKGLONGNAME} ${NIGHT}
rm -r /var/spool/pkg/${PKGNAME}

echo "     ${PKGLONGNAME} created" >> ${NIGHTLOG}
echo `date "+%d.%m.%y %T"`"   ${NAME} solaris package creation succeeded" >> ${NIGHTLOG}

# create documentation ?

cd ~
rm -rf ${BUILDDIR}
rm -rf ${INSTALLDIR}
