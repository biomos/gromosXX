#!/bin/sh
## Original script by Dug Song <dugsong@UMICH.EDU>, used for a chrooted
## postfix environment, adapted by Jasper Aukes <J.Aukes@azu.nl> to be used for
## chrooted package creation.
#
# Location: /root/bin/setup_chroot
#
# Usage:
#
# First, create a chrooted dir with needed files (some may be left out (tcsh
# etc), some might be missing on your system, some might be located elsewhere.
# This script is just a dirty hack and by NO means intelligent. I should
# improve it and write it in Perl when i have time.
#
# /root/bin/setup_chroot /tmp/PAKNAME
# 
# Then, tar zxvf your package into: /tmp/PAKNAME/tmp/
#
# cd /tmp/PAKNAME/tmp, configure and compile your package 
# DO NOT 'make install' just yet...
#
# Now, cd to /tmp and run:
#
# chroot PAKNAME /bin/sh        # Or /bin/tcsh if you can and if you prefer it
#
# Check if you're really in the chrooted environment (f.e. cat /etc/passwd, it
# shouldn't be there) :-)
#
# Now cd to your /tmp/package-version##         # fill in the right name
# and run a make install
#
# Exit your chrooted shell and start building the SUN package from the
# directory /tmp/PAKNAME/usr/local/     # Or /usr, or even /
#
# A nice page to see how this stage could see a happy ending is to be found at
# http://sunfreeware.com/pkgadd.html by Steven M. Christensen

PATH=/usr/bin:/sbin:/usr/sbin

# Create chroot'd area under Solaris 2.5.1 for postfix.
#
# Dug Song <dugsong@UMICH.EDU>

if [ $# -ne 1 ]; then
  echo "Usage: `basename $0` <directory>, e.g.: `basename $0` /tmp/pkgname" ; exit 1
fi

CHROOT=$1
  
# If CHROOT does not exist but parent does, create CHROOT
if [ ! -d ${CHROOT} ]; then
  # lack of -p below is intentional
  mkdir ${CHROOT}
fi
if [ ! -d ${CHROOT} -o "${CHROOT}" = "/" -o "${CHROOT}" = "/usr" ]; then
  echo "$0: bad chroot directory ${CHROOT}"
  exit 2
fi
for dir in etc/default etc/inet dev usr/bin lib usr/lib usr/share/lib/zoneinfo \
    usr/local net \
    tmp ; do
  if [ ! -d ${CHROOT}/${dir} ]; then mkdir -p ${CHROOT}/${dir} ; fi
done
ln -s usr/bin ${CHROOT}/bin
ln -s usr/bin ${CHROOT}/net/bin

# Set the right permissions
chmod -R 755 ${CHROOT}

# Copy some terminfo files
for term in v x ; do
  if [ ! -d ${CHROOT}/usr/share/lib/terminfo/${term} ]; then \
    mkdir -p ${CHROOT}/usr/share/lib/terminfo/${term} ; fi
  cp /usr/share/lib/terminfo/${term}/* ${CHROOT}/usr/share/lib/terminfo/${term}
  chmod 644 ${CHROOT}/usr/share/lib/terminfo/${term}/*
done

# AFS support.
if [ "`echo $CHROOT | cut -c1-4`" = "/afs" ]; then
  echo '\tCreating memory resident /dev...'
  mount -F tmpfs -o size=10 swap ${CHROOT}/dev
fi

# Setup /etc files.
cp /etc/nsswitch.conf ${CHROOT}/etc
cp /etc/netconfig /etc/resolv.conf ${CHROOT}/etc
cp /etc/default/init ${CHROOT}/etc/default
cp /etc/inet/services ${CHROOT}/etc/inet/services
ln -s /etc/inet/services ${CHROOT}/etc/services
cp /usr/share/lib/termcap ${CHROOT}/usr/share/lib
ln -s ${CHROOT}/usr/share/lib/termcap ${CHROOT}/etc/termcap
find ${CHROOT}/etc -type f -exec chmod 444 {} \;

# Most of the following are needed for basic operation, except
# for libnsl.so, nss_nis.so, libsocket.so, and straddr.so which are
# needed to resolve NIS names.
cp /usr/lib/ld.so.1 ${CHROOT}/usr/lib
cp /lib/ld.so.1 ${CHROOT}/lib

for lib in libc libdl libintl libmp libnsl libsocket libw libkstat \
    libcurses libkvm libelf libgen nss_nis nss_nisplus nss_dns nss_files; do
  cp /usr/lib/${lib}.so.1 ${CHROOT}/usr/lib
  rm -f ${CHROOT}/usr/lib/${lib}.so
  ln -s ./${lib}.so.1 ${CHROOT}/usr/lib/${lib}.so
done

for lib in straddr libmp; do
  cp /usr/lib/${lib}.so.2 ${CHROOT}/usr/lib
  rm -f ${CHROOT}/usr/lib/${lib}.so
  ln -s ./${lib}.so.2 ${CHROOT}/usr/lib/${lib}.so
done

chmod 555 ${CHROOT}/usr/lib/*

# Copy timezone database.
(cd ${CHROOT}/usr/share/lib/zoneinfo
  (cd /usr/share/lib/zoneinfo; find . -print | cpio -o) | cpio -imdu
  find . -print | xargs chmod 555
)

# Make device nodes. We need ticotsord, ticlts and udp to resolve NIS names.
for device in zero tcp udp ticotsord ticlts; do
  line=`ls -lL /dev/${device} | sed -e 's/,//'`
  major=`echo $line | awk '{print $5}'`
  minor=`echo $line | awk '{print $6}'`
  rm -f ${CHROOT}/dev/${device}
  mknod ${CHROOT}/dev/${device} c ${major} ${minor}
done
chmod 666 ${CHROOT}/dev/*

# Now copy some usefull binaries
for bin in expr ls dirname cp chmod rm mv sed mkdir grep find cat \
  true basename ln chown false cmp chgrp ; do
  cp /usr/bin/${bin} ${CHROOT}/usr/bin
done

cp /usr/ccs/bin/strip ${CHROOT}/usr/bin

cp /bin/sh ${CHROOT}/usr/bin

#for bin in install tar tcsh make ; do
#  cp /usr/local/bin/${bin} ${CHROOT}/usr/bin
#done

chmod 755 ${CHROOT}/usr/bin/*

exit 0


