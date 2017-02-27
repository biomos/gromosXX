#!/bin/sh

echo "preparing local settings"
echo ""

mkdir -p config
mv INSTALL INSTALL.bak
aclocal &&
libtoolize --copy --force&&
autoconf --force &&
autoheader --force &&
automake --add-missing --copy --force &&
autoheader --force ||
echo "setup failed. try doing it manually"

mv INSTALL.bak INSTALL

echo ""
echo "configure next"
echo ""
