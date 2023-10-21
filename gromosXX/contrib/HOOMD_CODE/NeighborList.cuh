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

// $Id: NeighborList.cuh 1983 2009-07-24 18:04:33Z akohlmey $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/cuda/NeighborList.cuh $
// Maintainer: joaander

#ifndef _NEIGHBORLIST_CUH_
#define _NEIGHBORLIST_CUH_

#include <stdio.h>
#include <cuda_runtime.h>

#include "ParticleData.cuh"

/*! \file NeighborList.cuh
	\brief Declares data structures and methods used by NeighborList and descendants
*/

/*! \ingroup data_structs
	@{
*/

/*! \defgroup gpu_data_structs GPU data structures
	\brief Data structures stored on the GPU
	\details See \ref page_dev_info for more information
*/

/*! @}
*/

//! Structure of arrays of the neighborlist data as it resides on the GPU
/*! The list is arranged as height rows of \a N columns with the given pitch.
	Each column i is the neighborlist for particle i. 
	\a n_neigh lists the number of neighbors in each column
	Each column has a fixed height \c height.
	
	\ingroup gpu_data_structs
*/
struct gpu_nlist_array
	{
	unsigned int *n_neigh;	//!< n_neigh[i] is the number of neighbors of particle with index i
	unsigned int *list;		//!< list[i*pitch + j] is the index of the j'th neighbor of particle i
	unsigned int height;	//!< Maximum number of neighbors that can be stored for any particle
	unsigned int pitch;		//!< width of the list in elements
	float4 *last_updated_pos;	//!< Holds the positions of the particles as they were when the list was last updated
	int *needs_update;		//!< Flag set to 1 when the neighbor list needs to be updated
	int *overflow;			//!< Flag set to 1 when the neighbor list overflows and needs to be expanded
	
	uint4 *exclusions;		//!< exclusions[i] lists all the particles that are to be excluded from being neighbors with [i] by index
#if defined(LARGE_EXCLUSION_LIST)
	uint4 *exclusions2;		//!< exclusions2[i] lists more the particles that are to be excluded from being neighbors with [i] by index
	uint4 *exclusions3;		//!< exclusions3[i] lists more the particles that are to be excluded from being neighbors with [i] by index
	uint4 *exclusions4;		//!< exclusions4[i] lists more the particles that are to be excluded from being neighbors with [i] by index
#endif
	};

//! Check if the neighborlist needs updating
cudaError_t gpu_nlist_needs_update_check(gpu_pdata_arrays *pdata, gpu_boxsize *box, gpu_nlist_array *nlist, float r_buffsq, int *result);

#endif
