#!/bin/sh

echo "preparing local settings"
echo ""

mkdir -p config
aclocal &&
libtoolize &&
autoconf &&
automake --add-missing ||
echo "setup failed. try doing it manually"

echo ""
echo "configure next"
echo ""
