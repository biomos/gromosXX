
#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "util/template_split.h"
#include "math/periodicity.h"

namespace interaction {
  // Specialize Impl for gpuBackend here:
  template<>
  class CUDA_Pairlist_Algorithm_Impl<util::gpuBackend> {
  public:
    CUDA_Pairlist_Algorithm_Impl() {
      DEBUG(0, "CUDA_Pairlist_Algorithm_Impl<util::gpuBackend> constructor");
    };
    
    int init(topology::Topology &topo, 
      configuration::Configuration &conf,
      simulation::Simulation &sim,
      std::ostream &os = std::cout,
      bool quiet = false) {
        // GPU-specific code
        DEBUG(0, "CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::init");
        return 0;
    };


    void prepare_cog(configuration::Configuration & conf,
                      topology::Topology & topo) {
        SPLIT_BOUNDARY(_prepare_cog, conf, topo);
    }
    

    /**
     * put the chargegroups into the box
     */
    template<math::boundary_enum b>
    void _prepare_cog(configuration::Configuration & conf,
                      topology::Topology & topo) {
        DEBUG(10, "putting chargegroups into box");
        math::Periodicity<b> periodicity(conf.current().box);
        periodicity.put_chargegroups_into_box(conf, topo);
    }

    math::CuFVArray m_cg_cog;
    float other_value;
    int value;
  };
}