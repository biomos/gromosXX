/**
 * @file cuda_kernel.h
 * contains Cuda_Kernel controlling the GPUs
 * contains function declarations used in GROMOSXX MD++
 */

#ifndef INCLUDED_CUDA_KERNEL_H
#define INCLUDED_CUDA_KERNEL_H

#ifdef HAVE_LIBCUDART
#include <cuda_runtime.h>
//#include <device_launch_parameters.h>
#include "gpu_status.h"
#else
#define gpu_status void
#endif

#include "io/message.h"

namespace topology { 
  class Topology;
}

namespace configuration { 
  class Configuration;
}

namespace simulation { 
  class Simulation;
}

namespace cudakernel {
    extern "C" {
        /**
         * @class CUDA_Kernel
         * performs operations with CUDA GPUs
         */
        class CUDA_Kernel {
        public:
            /**
             * Constructor
             */
            explicit CUDA_Kernel();

            /**
             * Destructor
             */
            ~CUDA_Kernel();

            void init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim);

            /**
             * Copy simulation parameters to GPU constant memory
             */
            void copy_parameters();
            
            /**
             * Estimate the pairlist size
             */
            void estimate_pairlist();
            
            void disabled() {
                io::messages.add("Compilation without CUDA support.",
                                "CUDA_Kernel", io::message::critical);
            };
        private:
            cudakernel::simulation_parameter param;
            #ifdef HAVE_LIBCUDART
            cudaDeviceProp device_properties;
            #endif
            topology::Topology * mytopo;
            configuration::Configuration * myconf;
            simulation::Simulation * mysim;
            int device_number;
            unsigned flags;
        };
    }
}
#endif
