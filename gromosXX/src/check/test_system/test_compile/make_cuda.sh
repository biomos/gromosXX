NAME="testbuild_cuda"

. ../gromos_settings.sh


cd ${GROMSRC}

./Config.sh


mkdir ${NAME}
cd ${NAME} 

../configure --with-cuda=/usr/local/cuda-11.4/  || $(echo "OHOH! ->conf" && exit 1);

make -j8  || $(echo "OHOH! -> MAKE" && exit 1);
make install  || $(echo "OHOH! -> MAKE install" && exit 1);

