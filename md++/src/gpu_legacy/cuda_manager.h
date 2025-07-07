/**
 * @file cuda_manager.h
 * contains CudaManager controlling the GPUs
 * contains function declarations used in GROMOSXX MD++
 */

#ifndef INCLUDED_CUDA_MANAGER_H
#define INCLUDED_CUDA_MANAGER_H

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

namespace cuda {
    /**
     * @class CudaManager
     * performs operations with CUDA GPUs
     */
    class CudaManager {
    public:
        /**
         * Get the instance object
         */
#ifdef HAVE_LIBCUDART
        CudaManager(std::shared_ptr<topology::Topology> topo,
                    std::shared_ptr<configuration::Configuration> conf,
                    std::shared_ptr<simulation::Simulation> sim);
#else
        CudaManager(std::shared_ptr<topology::Topology> topo,
                    std::shared_ptr<configuration::Configuration> conf,
                    std::shared_ptr<simulation::Simulation> sim) {
            io::messages.add(
                    "This binary was compiled without CUDA support, "
                    "but CUDA-dependent functionality was requested at runtime.",
                    "CudaManager", io::message::critical);
        }
#endif
        /**
         * disable copy
         */
        CudaManager(const CudaManager &) = delete;
        /**
         * disable assignment
         */
        CudaManager& operator=(const CudaManager &) = delete;

        /**
         * Destructor
         */
        ~CudaManager();

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
                            "CudaManager", io::message::critical);
        };
    protected:
    private:

        
        /**
         * CUDA-specific simulation parameters
         */
        cukernel::simulation_parameter param;

#ifdef HAVE_LIBCUDART
        cudaDeviceProp device_properties;
#endif
        /**
         * topology
         */
        std::shared_ptr<topology::Topology> mytopo;

        /**
         * configuration
         */
        std::shared_ptr<topology::Topology> myconf;

        /**
         * simulation
         */
        std::shared_ptr<topology::Topology> mysim;

        /**
         * Nonbonded_Parameters
         */
        //interaction::Nonbonded_Parameter * myparam;

        int device_number; // maybe we store the device numbers here and run all the GPUs from here
        unsigned flags;
    };
}
#endif
