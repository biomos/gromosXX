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

// $Id: ParticleData.cuh 1826 2009-04-27 22:18:31Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/cuda/ParticleData.cuh $
// Maintainer: joaander

#ifndef _PARTICLEDATA_CUH_
#define _PARTICLEDATA_CUH_

#include <stdio.h>
#include <cuda_runtime.h>

/*! \file ParticleData.cuh
 	\brief Declares GPU kernel code and data structure functions used by ParticleData
*/

//! Structure of arrays of the particle data as it resides on the GPU
/*! Stores pointers to the particles positions, velocities, acceleartions, and particle tags.
	Particle type information is most likely needed along with the position, so the type
	is encoded in the 4th float in the position float4 as an integer. Device code
	can decode this type data with __float_as_int();
	
	All the pointers in this structure are allocated on the device.
	
	This structure is about to be rewritten. Consider it being documented as poorly documented
	for now.
	
	\ingroup gpu_data_structs
*/
struct gpu_pdata_arrays
	{
	unsigned int N;	//!< Number of particles in the arrays
	unsigned int local_beg;	//!< Index of the first particle local to this GPU
	unsigned int local_num;	//!< Number of particles local to this GPU	
	
	float4 *pos;	//!< Particle position in \c x,\c y,\c z, particle type as an int in \c w
	float4 *vel;	//!< Particle velocity in \c x, \c y, \c z, nothing in \c w
	float4 *accel;	//!< Particle acceleration in \c x, \c y, \c z, nothing in \c w
	float *charge;	//!< Particle charge
	float *mass;	//!< Particle mass
	float *diameter;	//!< Particle diameter
	int4 *image;	//!< Particle box image location in \c x, c y, and \c z. Nothing in \c w.
	
	unsigned int *tag;	//!< Particle tag
	unsigned int *rtag;	//!< Particle rtag 
	};

//! Store the box size on the GPU
/*!	\note For performance reasons, the GPU code is allowed to assume that the box goes
	from -L/2 to L/2, and the box dimensions in this structure must reflect that.
	
	\ingroup gpu_data_structs
*/
struct gpu_boxsize
	{
	float Lx;	//!< Length of the box in the x-direction
	float Ly;	//!< Length of the box in the y-direction
	float Lz;	//!< Length of teh box in the z-direction
	float Lxinv;//!< 1.0f/Lx
	float Lyinv;//!< 1.0f/Ly
	float Lzinv;//!< 1.0f/Lz
	};
	
//! Helper kernel for un-interleaving data
cudaError_t gpu_uninterleave_float4(float *d_out, float4 *d_in, int N, int pitch);
//! Helper kernel for interleaving data
cudaError_t gpu_interleave_float4(float4 *d_out, float *d_in, int N, int pitch);

//! Generate a test pattern in the data on the GPU (for unit testing)
cudaError_t gpu_generate_pdata_test(const gpu_pdata_arrays &pdata);
//! Read from the pos texture and write into the vel array (for unit testing)
cudaError_t gpu_pdata_texread_test(const gpu_pdata_arrays &pdata); 

#endif	
