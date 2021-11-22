
. ../gromos_settings.sh

for x in $( ls ${GROMSRC}/testbuild* -d);
do 

	if [ ${x} != "*mpi*" ]; 
	then
		echo ${x};
		binDIR="${x}/bin"
		dirP=$(basename $x)
		dirP=${dirP/build_}

		echo ${dirP}
		rm -r ${dirP}
		cp -r template $dirP

		cd ${dirP}
		./job.sh ${binDIR}
		cd ..
        fi
done

