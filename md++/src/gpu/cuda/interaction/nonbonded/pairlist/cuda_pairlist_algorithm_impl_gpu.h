
#pragma once

namespace interaction {
  // Specialize Impl for gpuBackend here:
  template<>
  class CUDA_Pairlist_Algorithm_Impl<util::gpuBackend> {
  public:
    CUDA_Pairlist_Algorithm_Impl();
    
    int init(topology::Topology &topo, 
      configuration::Configuration &conf,
      simulation::Simulation &sim,
      std::ostream &os = std::cout,
      bool quiet = false);


    void prepare_cog(configuration::Configuration & conf,
                      topology::Topology & topo);

    /**
     * put the chargegroups into the box
     */
    template<math::boundary_enum b>
    void _prepare_cog(configuration::Configuration & conf,
                      topology::Topology & topo);

    math::CuVArray m_cg_cog;
    float other_value;
    int value;
  };
}