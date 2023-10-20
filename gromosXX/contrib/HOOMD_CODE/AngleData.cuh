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

// $Id: AngleData.cuh 1862 2009-06-08 21:07:24Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/cuda/AngleData.cuh $
// Maintainer: dnlebard

#ifndef _ANGLEDATA_CUH_
#define _ANGLEDATA_CUH_

#include <stdio.h>
#include <cuda_runtime.h>

/*! \file AngleData.cuh
 	\brief GPU data structures used in AngleData
*/

//! Angle data stored on the GPU
/*! gpu_angletable_array stores all of the angles between particles on the GPU.
	It is structured similar to gpu_nlist_array in that a single column in the list
	stores all of the angles for the particle associated with that column. 
	
	To access angle \em a of particle with local index \em i, use the following indexing scheme
	\code
	uint2 angle = angletable.angles[a*angletable.pitch + i]
	\endcode
	The particle with \b index (not tag) \em i is angled with particles \em angle.x
        and \em angle.y	with angle type \em angle.z. Each particle may have a different number of angles as
	indicated in \em n_angles[i].
	
	Only \a num_local angles are stored on each GPU for the local particles
	
	\ingroup gpu_data_structs
*/


struct gpu_angletable_array
	{
	unsigned int *n_angles;	//!< Number of angles for each particle
	uint4 *angles;			//!< angle list
	unsigned int height;	//!< height of the angle list
	unsigned int pitch;		//!< width (in elements) of the angle list

	//! Allocates memory
	cudaError_t allocate(unsigned int num_local, unsigned int alloc_height);
	
	//! Frees memory
	cudaError_t deallocate();
	};

#endif
