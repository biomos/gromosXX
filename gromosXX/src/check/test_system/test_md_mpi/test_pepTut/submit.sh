#Load Gromos Settings
. ../gromos_settings.sh

for x in $( ls ${GROMDIR}/testbuild* -d);
do 
	if [[ ${x} == *"testbuild"*"mpi"* ]];
	then
		echo "TEST BUILD: ${x}";

		binDIR="${x}/bin"
		dirP=$(basename $x)
		dirP=${dirP/build}

		echo "dir: ${dirP}"
		rm -rf ${dirP} #DELETE previous data
		cp -r templateMdMpi ${dirP}

		cd ${dirP}
		echo "./submit_all.sh ${binDIR}"
		./submit_all.sh ${binDIR} || echo "Failed: ${x} " >> Failed.log
		cd ..
	fi;
done

