for x in $(ls make*)
do
	./${x} || echo "FAILED IN : ${x}" && exit 1

done
