/**
 * @file interaction.h
 * interaction computation
 */
#ifndef CUKERNEL_INTERACTION
namespace cudakernel {
/**
  * calculate the forces, energies and virial and stores them in per atom arrays
  * @param[in] pl the pairlist used for the interaction computation
  * @param[in] dev_pos the positions
  * @param[in] dev_params the simulation parameters
  * @param[in] dev_lj_crf_params nonbonded parameters
  * @param[out] dev_for the forces
  * @param[out] dev_virial the 3x3 virial tensor for every atom
  * @param[out] dev_energy the LJ and CRF energies for every atom
  */
__global__ void kernel_CalcForces_Solvent(
        pairlist pl,
        float3 * dev_pos,
        cudakernel::simulation_parameter * dev_params,
        lj_crf_parameter * dev_lj_crf_params,
        float3 * dev_for,
        float9 * dev_virial,
        float2 * dev_energy,
        unsigned int num_of_gpus,
        unsigned int gpu_id);
}
#endif

