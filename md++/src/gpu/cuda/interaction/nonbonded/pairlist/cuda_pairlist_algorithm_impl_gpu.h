
#pragma once

namespace interaction {
  // Specialize Impl for gpuBackend here:
  template<>
  class CUDA_Pairlist_Algorithm_Impl<util::gpuBackend> : public CUDA_Pairlist_Algorithm_ImplBase {
    public:
      CUDA_Pairlist_Algorithm_Impl();
      
      int init(topology::Topology &topo, 
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os = std::cout,
        bool quiet = false);


      void prepare_cog(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

      /**
       * put the chargegroups into the box
       */
      template<math::boundary_enum b>
      void _prepare_cog(configuration::Configuration & conf,
                        topology::Topology & topo);
              

      void set_cutoff(double const cutoff_short, double const cutoff_long) override;

      /**
       * order atoms / chargegroups based on their grid cell?
       */
      void reorder(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

    private:
      /**
       * chargegroup center of geometry array.
       */
      math::CuVArray m_cg_cog;
      /**
       * chargegroup cell indices array.
       */
      gpu::cuvector<ushort4> m_cg_cells;
  };
}