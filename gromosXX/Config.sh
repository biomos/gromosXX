#!/bin/sh

echo "preparing local settings"
echo ""

mkdir -p config
aclocal &&
libtoolize --copy &&
autoconf &&
autoheader &&
automake --add-missing --copy &&
autoheader ||
echo "setup failed. try doing it manually"

echo ""
echo "configure next"
echo ""
