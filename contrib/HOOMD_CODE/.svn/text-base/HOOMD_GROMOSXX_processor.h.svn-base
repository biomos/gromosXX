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
