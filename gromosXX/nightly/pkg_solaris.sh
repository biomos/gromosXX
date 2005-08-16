#!/usr/bin/bash

echo "     preparing directory"	
mkdir -p /tmp/nightly_gromosXX_${BUILD}
cp ~/RELEASE/nightly/gromosXX-${VERSION}-${BUILD}.tar.gz /tmp/nightly_gromosXX_${BUILD}
cd /tmp/nightly_gromosXX_${BUILD}

${TAR} zxvf gromosXX-${VERSION}-${BUILD}.tar.gz > ${NIGHT}/log/pkg_solaris.log 2>&1
cd gromosXX-${VERSION}

echo "     configuring"
./configure ${CONFIGURE_REL_OPT} >> ${NIGHT}/log/pkg_solaris.log 2>&1 &&
    echo "     make" &&
    make >> ${NIGHT}/log/pkg_solaris.log 2>&1 &&
    echo "     enter fakeroot" &&
    fakeroot /usr/bin/bash nightly/pkg_solaris_fakeroot.sh

cd /tmp/gromosXX_pkg_solaris

cp /tmp/nightly_gromosXX_${BUILD}/gromosXX-${VERSION}/packaging/solaris.depend depend
cp /tmp/nightly_gromosXX_${BUILD}/gromosXX-${VERSION}/copyright .

ARCH=`uname -p`
DATE=`date +%d.%m.%y`

echo 'PKG="IGCgromosxx"' > pkginfo
echo "NAME=\"IGC GromosXX ${VERSION} ${BUILD}\"" >> pkginfo
echo "VERSION=\"${VERSION} ${BUILD}\"" >> pkginfo
echo "ARCH=\"${ARCH}\"" >> pkginfo
echo "CLASSES=\"none\"" >> pkginfo
echo "CATEGORY=\"develop\"" >> pkginfo
echo "VENDOR=\"IGC\"" >> pkginfo
echo "PSTAMP=\"${DATE}\"" >> pkginfo
echo "BASEDIR=/" >> pkginfo

pkgmk -r /tmp/gromosXX_pkg_solaris
pkgtrans -s /var/spool/pkg /tmp/IGCGromosXX-${VERSION}-${BUILD}-${OS}.pkg IGCgromosxx

mv /tmp/IGCGromosXX-${VERSION}-${BUILD}-${OS}.pkg ~/RELEASE/nightly/
rm -r /var/spool/pkg/IGCgromosxx

echo `date "+%d.%m.%y %T"`"     solaris package creation succeeded" >> ${NIGHTLOG}

# create documentation ?

cd ~
rm -rf /tmp/nightly_gromosXX_${BUILD} > /dev/null 2>&1
rm -rf /tmp/gromosXX_pkg_solaris

