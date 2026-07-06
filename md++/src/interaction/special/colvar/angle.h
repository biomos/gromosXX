#ifndef INCLUDED_ANGLE_COLVAR
#define INCLUDED_ANGLE_COLVAR

#include "colvar.h"
#include "../../../topology/topology.h"
#include "../../../util/virtual_atom.h"

namespace interaction {

class Angle_Colvar : public Colvar {
public:
  Angle_Colvar() : Colvar("Angle") {}

  topology::angle_restraint_struct *params;

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

private:
  std::vector<util::Virtual_Atom> m_atoms;
};

}

#endif
