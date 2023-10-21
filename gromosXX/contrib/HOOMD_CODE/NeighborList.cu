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

// $Id: NeighborList.cu 1826 2009-04-27 22:18:31Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/cuda/NeighborList.cu $
// Maintainer: joaander

#include "NeighborList.cuh"
#include "gpu_settings.h"

#ifdef WIN32
#include <cassert>
#else
#include <assert.h>
#endif

/*! \file NeighborList.cu
	\brief Defines data structures and methods used by NeighborList and descendants
*/

//! Checks if the neighbor list needs updating
/*! \param pdata Particle data to check
	\param nlist Current neighbor list build from that particle data
	\param r_buffsq A precalculated copy of r_buff*r_buff
	\param box Box dimensions for periodic boundary handling
	
	If any particle has moved a distance larger than r_buffsq since the last neighbor list update, 
	nlist.needs_update is set to 1.
*/
__global__ void gpu_nlist_needs_update_check_kernel(gpu_pdata_arrays pdata, gpu_nlist_array nlist, float r_buffsq, gpu_boxsize box)
	{
	// each thread will compare vs it's old position to see if the list needs updating
	// if that is true, write a 1 to nlist_needs_updating
	// it is possible that writes will collide, but at least one will succeed and that is all that matters
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int pidx = idx + pdata.local_beg;

	if (idx < pdata.local_num)
		{
		float4 cur_pos = pdata.pos[pidx];
		float4 last_pos = nlist.last_updated_pos[pidx];
		float dx = cur_pos.x - last_pos.x;
		float dy = cur_pos.y - last_pos.y;
		float dz = cur_pos.z - last_pos.z;
	
		dx = dx - box.Lx * rintf(dx * box.Lxinv);
		dy = dy - box.Ly * rintf(dy * box.Lyinv);
		dz = dz - box.Lz * rintf(dz * box.Lzinv);
	
		float drsq = dx*dx + dy*dy + dz*dz;

		if (drsq >= r_buffsq && pidx < pdata.N)
			{
			*nlist.needs_update = 1;
			}
		}
	}

//! Check if the neighborlist needs updating
/*! \param pdata Particle data to check
	\param box Box dimensions for periodic boundary handling
	\param nlist Current neighbor list build from that particle data
	\param r_buffsq A precalculated copy of r_buff*r_buff
	\param result Pointer to write the result to
	
	If any particle has moved a distance larger than r_buffsq since the last neighbor list update, 
	*result is set to 1. Otherwide *result is set to 0.
*/
cudaError_t gpu_nlist_needs_update_check(gpu_pdata_arrays *pdata, gpu_boxsize *box, gpu_nlist_array *nlist, float r_buffsq, int *result)
	{
	assert(pdata);
	assert(nlist);
	assert(result);
	
	// start by zeroing the value on the device
	*result = 0;
	cudaError_t error = cudaMemcpy(nlist->needs_update, result,
			sizeof(int), cudaMemcpyHostToDevice);
	
	// run the kernel
	int M = 256;
	dim3 grid( (pdata->local_num/M) + 1, 1, 1);
	dim3 threads(M, 1, 1);

	// run the kernel
	if (error == cudaSuccess)
		{
		gpu_nlist_needs_update_check_kernel<<< grid, threads >>>(*pdata, *nlist, r_buffsq, *box);
		if (!g_gpu_error_checking)
			{
			error = cudaSuccess;
			}
		else
			{
			cudaThreadSynchronize();
			error = cudaGetLastError();
			}
		}
	
	if (error == cudaSuccess)
		{
		error = cudaMemcpy(result, nlist->needs_update,
			sizeof(int), cudaMemcpyDeviceToHost);
		}
	return error;
	}

// vim:syntax=cpp
