
LD_LIBRARY_PATH=/opt/csw/lib
export LD_LIBRARY_PATH

make DESTDIR=/tmp/gromosXX_pkg_solaris install
cd /tmp/gromosXX_pkg_solaris
(echo 'i pkginfo'; echo 'i copyright'; echo 'i depend'; pkgproto .) | sed 's: none : none /:' > prototype

unset LD_LIBRARY_PATH
