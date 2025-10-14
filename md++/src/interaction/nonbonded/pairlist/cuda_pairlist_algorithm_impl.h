
#pragma once

#include "topology/topology.h"
#include "configuration/configuration.h"

#include "util/template_split.h"
#include "math/periodicity.h"

namespace interaction {
    class CUDA_Pairlist_Algorithm_ImplBase {
        public:
            virtual void set_cutoff(double const cutoff_short, double const cutoff_long) = 0;
        
            /**
             * squared shortrange cutoff.
             */
            double m_cutoff_short_2;

            /**
             * squared longrange cutoff.
             */
            double m_cutoff_long_2;

            /**
             * longrange cutoff.
             */
            double m_cutoff_long;

            /**
             * shortrange cutoff.
             */
            double m_cutoff_short;
    };

    template <typename Backend>
    class CUDA_Pairlist_Algorithm_Impl : public CUDA_Pairlist_Algorithm_ImplBase {
    public:
        CUDA_Pairlist_Algorithm_Impl() : CUDA_Pairlist_Algorithm_ImplBase() {};

        void set_cutoff(double const cutoff_short, double const cutoff_long) override
        {
            m_cutoff_long = cutoff_long;
            m_cutoff_short = cutoff_short;
            m_cutoff_short_2 = cutoff_short * cutoff_short;
            m_cutoff_long_2  = cutoff_long * cutoff_long;
        }

        void prepare_cog(configuration::Configuration & conf,
                        topology::Topology & topo,
                        simulation::Simulation & sim) {
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

        /**
         * order atoms / chargegroups based on their grid cell
         */
        void reorder(configuration::Configuration & conf,
                    topology::Topology & topo,
                    simulation::Simulation & sim) {
                        
                    };
    };
}

// force instantiation
// template class interaction::CUDA_Pairlist_Algorithm_Impl<util::cpuBackend>;

#ifdef USE_CUDA
  #include "gpu/cuda/interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl_gpu.h"
#endif
