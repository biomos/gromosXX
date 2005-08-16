#!/usr/bin/bash

echo "Debian package" >> ${NIGHT}/log/pkg_debian.log
echo "     preparing directory"	>> ${NIGHT}/log/pkg_debian.log
mkdir -p ${TMP}/nightly_gromosXX_${BUILD}
cp ~/RELEASE/nightly/gromosXX-${VERSION}-${BUILD}.tar.gz ${TMP}/nightly_gromosXX_${BUILD}
cd ${TMP}/nightly_gromosXX_${BUILD}

${TAR} zxvf gromosXX-${VERSION}-${BUILD}.tar.gz > ${NIGHT}/log/pkg_debian.log 2>&1
cd gromosXX-${VERSION}

export LD_LIBRARY_PATH=~/programs/x86/opt/intel/cc/9.0/lib

echo "     configuring" >> ${NIGHT}/log/pkg_debian.log
./configure ${CONFIGURE_REL_OPT} >> ${NIGHT}/log/pkg_debian.log 2>&1 &&
    echo "     make" >> ${NIGHT}/log/pkg_debian.log &&
    make >> ${NIGHT}/log/pkg_debian.log 2>&1 &&
    echo "     enter fakeroot" &&
    fakeroot /usr/bin/bash nightly/pkg_debian_fakeroot.sh

unset LD_LIBRARY_PATH

cd ${TMP}/gromosXX_pkg_debian

mkdir DEBIAN

ARCH=`uname -p`
DATE=`date +%d.%m.%y`

echo 'Package: IGCgromosxx' > DEBIAN/control
echo "Version: ${VERSION} ${BUILD}" >> DEBIAN/control
echo "Section: simulation" >> DEBIAN/control
echo "Priority: optional" >> DEBIAN/control
echo "Architecture: ${ARCH}" >> DEBIAN/control
echo "Depends: libc6" >> DEBIAN/control
echo "Depends: libgcc1" >> DEBIAN/control
echo "Depends: libstdc++5" >> DEBIAN/control
echo "Depends: libgsl0" >> DEBIAN/control
echo "Source: IGCgromosxx-src" >> DEBIAN/control
echo "Maintainer: Markus Christen <markus@igc.phys.chem.ethz.ch>" >> DEBIAN/control
echo "Description: GromosXX (Groningen Molecular Simulation)" >> DEBIAN/control

cd ..
dpkg --build gromosXX_pkg_debian

mv ${TMP}/IGCgromosXX.deb ~/RELEASE/nightly/IGCGromosXX-${VERSION}-${BUILD}-${OS}.deb
rm -r ${TMP}/IGCgromosxx.deb

echo `date "+%d.%m.%y %T"`"     debian package creation succeeded" >> ${NIGHTLOG}

# create documentation ?

cd ~
rm -rf ${TMP}/nightly_gromosXX_${BUILD} > /dev/null 2>&1
rm -rf ${TMP}/gromosXX_pkg_debian

