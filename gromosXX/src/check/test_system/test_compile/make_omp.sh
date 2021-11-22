NAME="testbuild_omp"

. ../gromos_settings.sh


cd ${GROMDIR}
./Config.sh

mkdir ${NAME}
cd ${NAME} 

../configure --enable-openmp  > ${NAME}.log  || $(echo "OHOH! ->conf" && exit 1);

make -j${compileCores}  > ${NAME}.log  || $(echo "OHOH! -> MAKE" && exit 1);
make install -j${compileCores}  > ${NAME}.log  || $(echo "OHOH! -> MAKE install" && exit 1);

