
. ../gromos_settings.sh

for x in $( ls ${GROMDIR}/testbuild* -d);
do 
	echo "${x}"
	if [[ ${x} == *"testbuild"*"mpi"* ]];
	then
		echo "TEST ${x}"
		binDIR="${x}/bin"
		dirP=$(basename $x)
		dirP=${dirP/build_}

		echo ${dirP}
		rm -r ${dirP}
		cp -r templateMdMpi $dirP

		cd ${dirP}
		./submit_all.sh ${binDIR}
		cd ..
	fi;
done

