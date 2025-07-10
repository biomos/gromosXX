/**
 * @file cuda_kernel.h
 * @author Poliak (peter.poliak@boku.ac.at)
 * @brief singleton controlling the GPUs
 * @version 0.1
 * @date 17.06.2023
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef INCLUDED_CUDA_KERNEL_H
#define INCLUDED_CUDA_KERNEL_H
/**
 * In all header files, #pragma once is used instead.
 * Although non-standard, these files are compilable
 * only using nvcc compiler, which is not about to
 * deprecate it.
 */

#ifdef HAVE_LIBCUDART
#include <cuda_runtime.h>
//#include <device_launch_parameters.h>
#include "gpu_status.h"
#include "lib/cuvector.h"
#include "lib/exclusions.h"

#else
#define gpu_status void
#endif

#include "io/message.h"

namespace math {
  class Box;
}

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

namespace cukernel {
    //extern "C" {
        /**
         * @class CUDA_Kernel
         * performs operations with CUDA GPUs
         */
        class CUDA_Kernel {
        public:
            /**
             * Get the instance object
             */
            static CUDA_Kernel* get_instance(const topology::Topology & topo,
                                             const configuration::Configuration & conf,
                                             const simulation::Simulation & sim);
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

            void init(const topology::Topology & topo,
                      const configuration::Configuration & conf,
                      const simulation::Simulation & sim);

            /**
             * Update nonbonded parameters
             */
            void update_nonbonded(interaction::Nonbonded_Parameter *np);

            /**
             * Copy parameters to constant device memory
             */
            void set_device_parameter() const;

            /**
             * Get parameters from device (for debugging)
             */
            cukernel::simulation_parameter get_device_parameter() const;

            /**
             * Copy simulation parameters to GPU constant memory
             */
            void copy_parameters();

            /**
             * Copy simulation box to GPU constant memory
             */
            int copy_box(math::Box & box);

            /**
             * Copy positions
             */
            int copy_positions(math::VArray & pos);
            
            /**
             * Estimate the pairlist size
             */
            void estimate_pairlist();
            
            /**
             * Calculate the nonbonded interactions
             */
            void calculate_interactions();

            /**
             * Update the pairlist
             */
            void update_pairlist(topology::Topology &topo,
                                 configuration::Configuration &conf,
                                 simulation::Simulation &sim);

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

            /**
             * The instance
             */
            static CUDA_Kernel * m_cuda_kernel;

            /**
             * Copy of simulation parameters transfered to device
             */
            cukernel::simulation_parameter param;

            /**
             * Copy of nonbonded exclusions transfered to device
             */
            cukernel::Exclusions exclusions;

#ifdef HAVE_LIBCUDART
            /**
             * Properties of the GPU
             */
            cudaDeviceProp device_prop;

            /**
             * positions for the GPU
             */
            cuvector<float3> m_pos;

            /**
             * forces for the GPU
             */
            cuvector<float3> m_force;
#endif
            /**
             * topology
             */
            const topology::Topology * mytopo;

            /**
             * configuration
             */
            const configuration::Configuration * myconf;

            /**
             * simulation
             */
            const simulation::Simulation * mysim;

            /**
             * Nonboded_Parameters
             */
            //interaction::Nonbonded_Parameter * myparam;

            int device_number; // maybe we store the device numbers here and run all the GPUs from here
            unsigned flags;
        };
    //}
}
#endif
