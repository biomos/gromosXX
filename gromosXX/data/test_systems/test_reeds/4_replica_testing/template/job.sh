
out_prefix="${1}"
gromosDIR=/cluster/work/igc/bschroed/gromos/gromos_compiled/gromosXX_more_Parallel_newEnvMPI2/bin

rm ${out_prefix}*
mpirun ${gromosDIR}/repex_mpi \
	@topo input/PNMT_9lig_water.top \
	@conf input/coord/REEDS_eoff_run.cnf \
	@input input/repex_eoff.imd \
	@pttopo input/PNMT_9lig_water.ptp \
	@distrest input/PNMT_9lig_water_disres.dat \
	@fin ${out_prefix}.cnf \
	@trc ${out_prefix}.trc \
	@tre ${out_prefix}.tre \
	@repout ${out_prefix}.dat \
	@repdat ${out_prefix}.dat >> ${out_prefix}.omd

