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

namespace interaction {
  class Nonbonded_Parameter;
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
             * Get the instance object
             */
            static CUDA_Kernel *get_instance(topology::Topology & topo,
                                             configuration::Configuration & conf,
                                             simulation::Simulation & sim);
            /**
             * disable copy
             */
            CUDA_Kernel(const CUDA_Kernel &) = delete;
            /**
             * disable assignment
             */
            void operator=(const CUDA_Kernel &) = delete;

            /**
             * Destructor
             */
            ~CUDA_Kernel();

            void init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim);

            /**
             * Update nonbonded parameters
             */
            void update_nonbonded(interaction::Nonbonded_Parameter *np);

            /**
             * Copy constants to device symbols
             */
            void sync_symbol();

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
        protected:
        private:
            /**
             * Constructor as singleton
             */
            CUDA_Kernel();

            cudakernel::simulation_parameter param;

#ifdef HAVE_LIBCUDART
            cudaDeviceProp device_properties;
#endif
            /**
             * The instance
             */
            static CUDA_Kernel * m_cuda_kernel;

            /**
             * topology
             */
            topology::Topology * mytopo;

            /**
             * configuration
             */
            configuration::Configuration * myconf;

            /**
             * simulation
             */
            simulation::Simulation * mysim;

            /**
             * Nonboded_Parameters
             */
            //interaction::Nonbonded_Parameter * myparam;

            int device_number; // maybe we store the device numbers here and run all the GPUs from here
            unsigned flags;
        };
    }
}
#endif
