#ifndef INCLUDED_COORDNUM_COLVAR
#define INCLUDED_COORDNUM_COLVAR

#include "colvar.h"
#include "../../../topology/topology.h"

namespace interaction {

class Coordnum_Colvar : public Colvar {
public:
  Coordnum_Colvar() : Colvar("CoordinationNumber") {}

  topology::coordnum_restraint_struct *params;

  int init(
    topology::Topology &topo,
    configuration::Configuration &conf,
    simulation::Simulation &sim,
    std::ostream &os,
    bool quiet);

  int calculate(
    topology::Topology &topo,
    configuration::Configuration &conf,
    simulation::Simulation &sim);
};

}

#endif
