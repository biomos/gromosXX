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
 * @file pairlist.h
 * pairlist computation
 */
#ifndef CUKERNEL_PAIRLIST
#define CUKERNEL_PAIRLIST
namespace cudakernel {
/**
 * @struct pairlist
 *  the holds the pairlist
 */
struct pairlist {
  /**
   * the elements in the pairlist [i*pitch + j]
   */
  unsigned int *list;
  /**
   * the number of neighbors in the list
   */
  unsigned int *num_neighbors;
  /**
   * the maximal number of neighbors that can be stored in the list
   */
  unsigned int max_size;
  /**
   * the pitch (from cudaMallocPitch)
   */
  unsigned int pitch;
  /**
   * bool overflow
   */
  bool *overflow;
};

#ifndef CUKERNEL_TYPES_ONLY
/**
 * free a pairlist object
 * @param pl the pairlist you want to free on the device
 */
void free_pairlist(pairlist &pl);
/**
 * allocate a pairlist object on the device
 * @param pl the pl object to hold the pointers
 * @param size the number of solvent molecules
 * @param max_neighbors the maximum number of neighbors the pairlist can hold
 */
void allocate_pairlist(pairlist &pl, unsigned int size, unsigned int max_neighbors);

/**
 * computes the pairlist on the GPU
 * @param[in] dev_param the simulation parameters
 * @param[in] dev_pos the positions
 * @param[out] pl_short the short-range pairlist
 * @param[out] pl_long the long-range pairlist
 */
__global__ void kernel_CalcPairlist(
        cudakernel::simulation_parameter * dev_params,
        float3 * dev_pos,
        pairlist pl_short,
        pairlist pl_long,
        unsigned int num_of_gpus,
        unsigned int gpu_id);
#endif
}
#endif

