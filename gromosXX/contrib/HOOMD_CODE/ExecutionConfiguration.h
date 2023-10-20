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

// $Id: ExecutionConfiguration.h 1919 2009-06-19 13:53:18Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/data_structures/ExecutionConfiguration.h $
// Maintainer: joaander

#ifndef __EXECUTION_CONFIGURATION__
#define __EXECUTION_CONFIGURATION__

#include "GPUWorker.h"

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

/*! \file ExecutionConfiguration.h
	\brief Declares ExecutionConfiguration and related classes
*/

//! Defines the execution configuration for the simulation
/*! \ingroup data_structs
	ExecutionConfiguration is a data structure needed to support HOOMD's unorthodox
	multi-GPU usage model. One GPUWorker thread is launched for each GPU that the run
	will execute on. The list of worker threads is maintained in ExecutionConfiguration::gpu.
	All calls to these GPUs must go through the corresponding worker thread in the list
	(see GPUWorker for examples).
	
	A few handy methods are defined to make working with more than one GPU simpler:
		- syncAll() calls sync() on each GPUWorker thread
		- callAll() makes a call() on each GPUWorker thread
		- tagAll()  tags the current position on each GPUWorker thread

	The execution configuration is determined at the beginning of the run and must 
	remain static for the entire run. It can be accessed from the ParticleData of the
	system. As to which particles/forces/etc go to which GPU in a multi-GPU environment,
	that is not determined here. See ParticleData and ForceCompute for a basic rundown
	on how that is broken up.	

	The execution mode is specified in exec_mode. This is only to be taken as a hint,
	different compute classes are free to execute on either the CPU or GPU. The only
	requirement is that those executing on the GPU must use the gpu workers in the vector.
*/
struct ExecutionConfiguration
	{
	//! Simple enum for the execution modes
	enum executionMode
		{
		GPU,	//!< Execute on the GPU
		CPU	//!< Execute on the CPU
		};
	
	//! Default constructor
	ExecutionConfiguration(bool min_cpu=false, bool ignore_display=false);
	
	//! Single GPU selection constructor
	ExecutionConfiguration(executionMode mode, bool min_cpu=false, bool ignore_display=false);
	
	//! Excplicit GPU selection constructor
	ExecutionConfiguration(const std::vector<int>& gpu_ids, bool min_cpu=false, bool ignore_display=false);
	
	executionMode exec_mode;	//!< Execution mode specified in the constructor
	
	#ifdef ENABLE_CUDA
	//! Sets tags for all GPUWorkers
	void tagAll(const std::string &file, unsigned int line) const;
	
	//! Syncs all GPUWorkers
	void syncAll() const;
	
	//! Calls a function on all GPUs
	void callAll(const boost::function< cudaError_t (void) > &func) const;
	
	std::vector< boost::shared_ptr<GPUWorker> > gpu;	//!< GPUs to execute on
	
	static int getDefaultGPU();	//!< returns the default GPU to run on
	static std::vector< int > getDefaultGPUList();	//!< returns the list of default GPUs to run on

	//! Checks all GPUs in the execution configuration to see if they meet the CUDA_ARCH min req.
	void checkCudaArch();
	
	//! Get the compute capability of the GPU that we are running on
	std::string getComputeCapability();
	
	private:
		//! Actually initializes the workers with the given list of GPUs
		void initializeGPUs(const std::vector<int>& gpu_ids, bool min_cpu);
		
		//! Print out stats on the chosen GPUs
		void printGPUStats();
		
		//! Scans through all GPUs reported by CUDA and marks if they are available
		void scanGPUs(bool ignore_display);
		
		//! Returns true if the given GPU is available for computation
		bool isGPUAvailable(int gpu_id);
		
		//! Returns the count of capable GPUs
		int getNumCapableGPUs();
		
		//! Return the number of GPUs that can be checked for availability
		unsigned int getNumTotalGPUs() { return m_gpu_available.size(); }
		
		std::vector< bool > m_gpu_available;    //!< true if the GPU is avaialble for computation, false if it is not
		bool m_system_compute_exclusive;        //!< true if every GPU in the system is marked compute-exclusive
		std::vector< int > m_gpu_list;          //!< A list of capable GPUs listed in priority order
	#endif
	};

	
#endif
