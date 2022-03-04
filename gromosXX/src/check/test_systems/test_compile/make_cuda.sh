NAME="testbuild_cuda"

. ../gromos_settings.sh


cd ${GROMDIR}
./Config.sh


mkdir ${NAME}
cd ${NAME} 


../configure --with-cuda=/usr/local/cuda-11.4/  || $(echo "OHOH! ->conf" && exit 1);

make -j${compileCores}  || $(echo "OHOH! -> MAKE" && exit 1);
make install -j${compileCores}   || $(echo "OHOH! -> MAKE install" && exit 1);

