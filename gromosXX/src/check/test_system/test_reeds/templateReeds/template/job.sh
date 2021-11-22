
out_prefix="${1}"
gromosBIN=${3}

rm -f ${out_prefix}*

mpirun -n ${2} ${gromosBIN}/repex_mpi \
	@topo input/PNMT_9lig_water.top \
	@conf input/coord/REEDS_eoff_run.cnf \
	@input input/repex_eoff.imd \
	@pttopo input/PNMT_9lig_water.ptp \
	@distrest input/PNMT_9lig_water_disres.dat \
	@fin ${out_prefix}.cnf \
	@trc ${out_prefix}.trc \
	@tre ${out_prefix}.tre \
	@repout ${out_prefix}.dat \
	@verb re:5 \
	@repdat ${out_prefix}.dat >> ${out_prefix}.omd

