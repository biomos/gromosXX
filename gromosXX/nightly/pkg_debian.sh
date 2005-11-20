#!/bin/bash

SCRIPT=pkg_debian
LOG=${NIGHT}/log/${NAME}_${SCRIPT}.log
BUILDDIR=${TMP}/nightly_${NAME}_${BUILD}
INSTALLDIR=${TMP}/${NAME}_${SCRIPT}

cat /dev/null > ${LOG}

echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation started" >> ${NIGHTLOG}

echo "creating a debian package" >> ${LOG}
echo "     preparing directory" >> ${LOG}
ok=1
export PATH=/usr/local/bin:/bin:/usr/bin:/sbin:/usr/sbin:/home/markus/programs/x86/bin:/home/markus/programs/x86/mpich2/bin:/home/markus/programs/x86/opt/intel/cc/9.0/bin:.

rm -rf ${BUILDDIR}
mkdir -p ${BUILDDIR} || ok=0
cp ${NIGHT}/${NAME}-${VERSION}-${BUILD}.tar.gz ${BUILDDIR} || ok=0
cd ${BUILDDIR} || ok=0

if [ ${ok} == 0 ] ; then
    echo "preparing directory failed"
    echo "preparing directory failed" >> ${LOG}
    echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation failed" >> ${NIGHTLOG}
    exit 1
fi

echo "    tar" >> ${LOG}
${TAR} zxvf ${NAME}-${VERSION}-${BUILD}.tar.gz > /dev/null 2>&1 || ok=0
cd ${NAME}-${VERSION} || ok=0

if [ ${ok} == 0 ] ; then
    echo "${TAR} zxvf ${NAME}-${VERSION}-${BUILD}.tar.gz failed"
    echo "${TAR} zxvf ${NAME}-${VERSION}-${BUILD}.tar.gz failed" >> ${LOG}
    echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation failed" >> ${NIGHTLOG}
    exit 1
fi

export LD_LIBRARY_PATH=~/programs/x86/opt/intel/cc/9.0/lib

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
echo "     fakeroot ${BASH} ${NIGHTHOME}/${SCRIPT}_fakeroot.sh" >> ${LOG}
fakeroot ${BASH} ${NIGHTHOME}/${SCRIPT}_fakeroot.sh || ok=0

if [ ${ok} == 0 ] ; then
    echo "enter fakeroot failed"
    echo "enter fakeroot failed" >> ${LOG}
    echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation failed" >> ${NIGHTLOG}
    exit 1
fi

cd ${INSTALLDIR} || ok=0
if [ ${ok} == 0 ] ; then
    echo "cd ${INSTALLDIR} failed"
    echo "cd ${INSTALLDIR} failed" >> ${LOG}
    echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation failed" >> ${NIGHTLOG}
    exit 1
fi

unset LD_LIBRARY_PATH

cd ${INSTALLDIR} || ok=0
if [ ${ok} == 0 ] ; then
    echo "cd ${INSTALLDIR} failed"
    echo "cd ${INSTALLDIR} failed" >> ${LOG}
    echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation failed" >> ${NIGHTLOG}
    exit 1
fi

mkdir DEBIAN || ok=0

ARCH=`uname -m`
DATE=`date +%d.%m.%y`

echo "Package: ${PKGNAME}" > DEBIAN/control || ok=0
echo "Version: ${VERSION}_${BUILD}" >> DEBIAN/control || ok=0
echo "Section: simulation" >> DEBIAN/control || ok=0
echo "Priority: optional" >> DEBIAN/control || ok=0
echo "Architecture: ${ARCH}" >> DEBIAN/control || ok=0
cat ${BUILDDIR}/${NAME}-${VERSION}/packaging/debian.depend >> DEBIAN/control || ok=0
echo "Source: ${PKGNAME}-src" >> DEBIAN/control || ok=0
echo "Maintainer: ${MAINTAINER}" >> DEBIAN/control || ok=0
echo "Description: ${DESCRIPTION}" >> DEBIAN/control || ok=0

if [ ${ok} == 0 ] ; then
    echo "creating control file failed"
    echo "creating control file failed" >> ${LOG}
    echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation failed" >> ${NIGHTLOG}
    exit 1
fi

cd ..
dpkg --build ${NAME}_${SCRIPT} >> ${LOG} 2>&1 || ok=0

if [ ${ok} == 0 ] ; then
    echo "dpkg --build ${NAME}_${SCRIPT} failed"
    echo "dpkg --build ${NAME}_${SCRIPT} failed" >> ${LOG}
    cat ${NAME}_${SCRIPT}/DEBIAN/control >> ${LOG}
    echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation failed" >> ${NIGHTLOG}
    exit 1
fi

mv ${TMP}/${NAME}_${SCRIPT}.deb ${NIGHT}/${PKGNAME}-${VERSION}-${BUILD}-${OS}.deb || ok=0
if [ ${ok} == 0 ] ; then
    echo "could not copy package"
    echo "could not copy package" >> ${LOG}
    echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation failed" >> ${NIGHTLOG}
    exit 1
fi

echo `date "+%d.%m.%y %T"`"   ${PKGNAME}-${VERSION}-${BUILD}-${OS}.deb created" >> ${NIGHTLOG}
echo `date "+%d.%m.%y %T"`"   ${NAME} debian package creation succeeded" >> ${NIGHTLOG}

# create documentation ?

cd ~

rm -rf ${BUILDDIR}
rm -rf ${INSTALLDIR}
