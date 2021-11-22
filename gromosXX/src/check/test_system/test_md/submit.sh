#Load Gromos Settings
. ../gromos_settings.sh

for x in $( ls ${GROMDIR}/testbuild* -d);
do 
		echo "TEST BUILD: ${x}";
		#path juggeling
		binDIR="${x}/bin"
		dirP=$(basename $x)
		dirP=${dirP/build}

		#build testFolder
		echo -e "\tgo to: ${dirP}"
		rm -rf ${dirP}	#remove possible old files
		cp -r template $dirP

		#execute test
		cd ${dirP}
		./job.sh ${binDIR} || echo "Failed with: ${binDIR}" >> FailedRund.out
		cd ..
done

