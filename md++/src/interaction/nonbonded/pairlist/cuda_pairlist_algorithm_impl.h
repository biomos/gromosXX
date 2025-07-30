
#pragma once

#include "topology/topology.h"
#include "configuration/configuration.h"

#include "util/template_split.h"
#include "math/periodicity.h"

namespace interaction {
    template <typename Backend>
    class CUDA_Pairlist_Algorithm_Impl {
    public:
        CUDA_Pairlist_Algorithm_Impl() {};

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
    };
}

// force instantiation
// template class interaction::CUDA_Pairlist_Algorithm_Impl<util::cpuBackend>;

#ifdef USE_CUDA
  #include "gpu/cuda/interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl_gpu.h"
#endif
