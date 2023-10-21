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

// $Id: BondData.cu 1826 2009-04-27 22:18:31Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/cuda/BondData.cu $
// Maintainer: joaander

#include "BondData.cuh"

#ifdef WIN32
#include <cassert>
#else
#include <assert.h>
#endif

/*! \file BondData.cu
 	\brief Defines the data structures for storing bonds on the GPU.
*/

/*! \pre no allocations have been performed or deallocate() has been called after a previous allocate()
	\post Memory for \a n_bonds and \a bonds is allocated on the device
	\param num_local Number of particles local to the GPU on which this is being called
	\param alloc_height Number of bonds to allocate for each particle
	\note allocate() \b must be called on the GPU it is to allocate data on
*/
cudaError_t gpu_bondtable_array::allocate(unsigned int num_local, unsigned int alloc_height)
	{
	// sanity checks
	assert(n_bonds == NULL);
	assert(bonds == NULL);
	
	height = alloc_height;
		
	// allocate n_bonds and check for errors
	cudaError_t error = cudaMalloc((void**)((void*)&n_bonds), num_local*sizeof(unsigned int));
	if (error != cudaSuccess)
		return error;
	
	error = cudaMemset((void*)n_bonds, 0, num_local*sizeof(unsigned int));
	if (error != cudaSuccess)
		return error;
	
	// cudaMallocPitch fails to work for coalesced reads here (dunno why), need to calculate pitch ourselves
	// round up to the nearest multiple of 32
	pitch = (num_local + (32 - num_local & 31));
	error = cudaMalloc((void**)((void*)&bonds), pitch * height * sizeof(uint2));
	if (error != cudaSuccess)
		return error;	
	
	error = cudaMemset((void*)bonds, 0, pitch * height * sizeof(uint2));
	if (error != cudaSuccess)
		return error;
		
	// all done, return success
	return cudaSuccess;	
	}
	
/*! \pre allocate() has been called
	\post Memory for \a n_bonds and \a bonds is freed on the device
	\note deallocate() \b must be called on the same GPU as allocate()
*/
cudaError_t gpu_bondtable_array::deallocate()
	{
	// sanity checks
	assert(n_bonds != NULL);
	assert(bonds != NULL);
	
	// free the memory
	cudaError_t error = cudaFree((void*)n_bonds);
	n_bonds = NULL;
	if (error != cudaSuccess)
		return error;
		
	error = cudaFree((void*)bonds);
	bonds = NULL;
	if (error != cudaSuccess)
		return error;
	
	// all done, return success
	return cudaSuccess;
	}
