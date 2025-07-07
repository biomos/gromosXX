/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file interaction.h
 * interaction computation
 */
#ifndef INCLUDED_CUKERNEL_INTERACTION_H
#define INCLUDED_CUKERNEL_INTERACTION_H

#define NUM_THREADS_PER_BLOCK_FORCES 192

#define NUM_THREADS_PER_BLOCK_REDUCTION_FLOAT2 512
#define NUM_ELEMENTS_PER_THREAD_FLOAT2 4
#define NUM_THREADS_PER_BLOCK_REDUCTION_FLOAT9 128
#define NUM_ELEMENTS_PER_THREAD_FLOAT9 4

namespace cukernel {
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
        cukernel::simulation_parameter * dev_params,
        lj_crf_parameter * dev_lj_crf_params,
        float3 * dev_for,
        float9 * dev_virial,
        float2 * dev_energy,
        unsigned int num_of_gpus,
        unsigned int gpu_id);
}
#endif

