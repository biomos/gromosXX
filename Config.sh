#!/bin/sh

echo "preparing local settings"
echo ""

mkdir -p config
aclocal &&
libtoolize --copy --force&&
autoconf --force &&
autoheader --force &&
automake --add-missing --copy --force &&
autoheader --force ||
echo "setup failed. try doing it manually"

svn revert INSTALL

echo ""
echo "configure next"
echo ""
