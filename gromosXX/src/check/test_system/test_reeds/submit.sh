
. ../../gromos_settings.sh

for x in $( ls ${GROMSRC}/testbuild* -d);
do 
	echo ${x};
	binDIR="${x}/bin"
	dirP=$(basename $x)
	dirP=${dirP/build_}

	echo ${dirP}
	rm -r ${dirP}
	cp -r templateReeds $dirP

	cd ${dirP}
	./submit_all.sh ${binDIR}
	cd ..

done

