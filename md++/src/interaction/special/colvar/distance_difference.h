/** Signed distance-difference collective variable: |A-B| - |C-D|. */
#ifndef INCLUDED_DISTANCE_DIFFERENCE_COLVAR_H
#define INCLUDED_DISTANCE_DIFFERENCE_COLVAR_H

#include "colvar.h"
#include "../../../topology/topology.h"

namespace interaction {

class Distance_Difference_Colvar : public Colvar {
public:
  Distance_Difference_Colvar() : Colvar("DistanceDifference"), params(NULL) {}

  int init(topology::Topology &topo, configuration::Configuration &conf,
           simulation::Simulation &sim, std::ostream &os = std::cout,
           bool quiet = false);
  int calculate(topology::Topology &topo, configuration::Configuration &conf,
                simulation::Simulation &sim);

  topology::distance_difference_restraint_struct *params;
};

} // namespace interaction
#endif
