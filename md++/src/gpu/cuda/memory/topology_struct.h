#pragma once

// #include "types.h"
// #include "cuvector.h"

namespace topology {
  class Topology;
}

namespace gpu {
  /**
   * @brief Holds GPU-side view of topology data.
   * 
   * All arrays are in single precision (float or float3) to match GPU performance requirements.
   * We assume constants, so we copy only once in the beginning
   */
  struct topology_struct {
    void*      memory_block; // base pointer for deallocation
    // Atom codes, masses, charges
    int*       iac;
    float*     mass;
    float*     inverse_mass;
    float*     charge;
    int*       chargegroup;
    unsigned   num_solute_chargegroups;
    unsigned   num_solute_molecules;
    unsigned   num_atoms;
    unsigned   num_chargegroups;

    topology_struct(topology::Topology& topo);
    ~topology_struct();

    void update(topology::Topology& topo);
  };

} // namespace configuration