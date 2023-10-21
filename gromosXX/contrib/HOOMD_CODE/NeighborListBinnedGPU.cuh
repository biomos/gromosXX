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

/*
Highly Optimized Object-Oriented Molecular Dynamics (HOOMD) Open
Source Software License
Copyright (c) 2008 Ames Laboratory Iowa State University
All rights reserved.

Redistribution and use of HOOMD, in source and binary forms, with or
without modification, are permitted, provided that the following
conditions are met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names HOOMD's
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND
CONTRIBUTORS ``AS IS''  AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS  BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Id: NeighborListBinnedGPU.cuh 1826 2009-04-27 22:18:31Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/cuda/NeighborListBinnedGPU.cuh $
// Maintainer: joaander

#ifndef _NEIGHBORLISTBINNEDGPU_CUH_
#define _NEIGHBORLISTBINNEDGPU_CUH_

#include <stdio.h>
#include <cuda_runtime.h>

#include "NeighborList.cuh"

/*! \file NeighborListBinnedGPU.cuh
	\brief Declares data structures and methods used by BinnedNeighborListGPU
*/

//! Structure of arrays storing the bins particles are placed in on the GPU
/*! This structure is in a current state of flux. Consider it documented as being
	poorly documented :)
	
	These are 4D arrays with indices i,j,k,n. i,j,k index the bin and each goes from 0 to Mx-1,My-1,Mz-1 respectively.
	Index into the data with idxdata[i*Nmax*Mz*My + j*Nmax*Mz + k*Nmax  + n] where n goes from 0 to Nmax - 1.
	
	\ingroup gpu_data_structs
*/
struct gpu_bin_array
        {
        unsigned int Mx;	//!< X-dimension of the cell grid
        unsigned int My;	//!< Y-dimension of the cell grid
        unsigned int Mz;	//!< Z-dimension of the cell grid
        unsigned int Nmax;	//!< Maximum number of particles each cell can hold
        unsigned int Nparticles;		//!< Total number of particles binned
        unsigned int coord_idxlist_width;	//!< Width of the coord_idxlist data
	
		unsigned int *bin_size;	//!< \a nbins length array listing the number of particles in each bin
        unsigned int *idxlist;	//!< \a Mx x \a My x \a Mz x \a Nmax 4D array holding the indices of the particles in each cell
		cudaArray *idxlist_array;	//!< An array memory copy of \a idxlist for 2D texturing
		
		float4 *coord_idxlist;	//!< \a Mx x \a My x \a Mz x \a Nmax 4D array holding the positions and indices of the particles in each cell (x,y,z are position and w holds the index)
		cudaArray *coord_idxlist_array;	//!< An array memory copy of \a coord_idxlist for 2D texturing
		
		cudaArray *bin_adj_array;	//!< bin_adj_array holds the neighboring bins of each bin in memory (x = idx (0:26), y = neighboring bin memory location)
		};
		

//! Take the idxlist and generate coord_idxlist
cudaError_t gpu_nlist_idxlist2coord(gpu_pdata_arrays *pdata, gpu_bin_array *bins, int curNmax, int block_size);

//! Kernel driver for GPU computation in BinnedNeighborListGPU
cudaError_t gpu_compute_nlist_binned(const gpu_nlist_array &nlist, const gpu_pdata_arrays &pdata, const gpu_boxsize &box, const gpu_bin_array &bins, float r_maxsq, int curNmax, int block_size, bool ulf_workaround);

#endif

