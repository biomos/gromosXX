for x in $(ls make*)
do
	./${x} || echo "FAILED IN : ${x}" > failedComp.out

done
