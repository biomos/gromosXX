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
 * GROMOSXX Interface to HOOMD Code
 * Wrapper for ExecutionConfiguration
 *
 * See documentation lyx/pdf
 */

#ifdef HAVE_HOOMD
#ifndef INCLUDED_HOOMD_GROMOSXX_PROCESSOR_H
#define INCLUDED_HOOMD_GROMOSXX_PROCESSOR_H

// Same as ones used in HOOMD_CODE/Makefile's CXXFLAGS. Change these also.
#define ENABLE_STATIC
#define SINGLE_PRECISION
#define CUDA_ARCH 10 // If the GPU architecture is CUDA Compute Capability 1.0, other options are 11, 12, 13
#define ENABLE_CUDA
#define LARGE_EXCLUSION_LIST

#include <ExecutionConfiguration.h>

namespace processor {
  enum which {
		 CPU,
		 GPUs
  };

  /**
   * @class Processor
   * defines which processors to use: CPU only or GPUs.
   */
  class Processor
  {
  public:
    ExecutionConfiguration *exec_conf;
	enum which chosen;
	
	/**
	 * Constructors
	 */
	Processor(enum which choice = GPUs) {
		switch (choice) {
			case CPU:
				exec_conf = new ExecutionConfiguration(ExecutionConfiguration::CPU);
				break;
			case GPUs:
				exec_conf = new ExecutionConfiguration(); 
				break;
		}
		chosen = choice;
	}

	/**
	 * Destructor
	 */
	~Processor() {
		delete exec_conf;
	}
  };
} // processor

#endif // INCLUDED_HOOMD_GROMOSXX_INTERFACE_H
#endif // HAVE_HOOMD
