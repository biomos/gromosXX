echo "Building"

mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release
make -j $(nproc)
#mkdir -p $PREFIX/bin
make install
