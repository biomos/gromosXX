#!/bin/sh

echo "preparing local settings"
echo ""

if [ -f src/interaction/nonbonded/split/Makefile.inc ] ; then
	echo "nonbonded split up ok"
	echo ""
else
	echo "nonbonded split up not ok"
	echo "see INSTALL about nonbonded split-up"
	echo ""
	exit
fi

mkdir -p config
aclocal &&
libtoolize --copy &&
autoconf &&
autoheader &&
automake --add-missing --copy --foreign ||
echo "setup failed. try doing it manually"

echo ""
echo "configure next"
echo ""
