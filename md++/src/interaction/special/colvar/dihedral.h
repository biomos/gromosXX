#ifndef INCLUDED_DIHEDRAL_COLVAR
#define INCLUDED_DIHEDRAL_COLVAR

#include "colvar.h"
#include "../../../topology/topology.h"
#include "../../../util/virtual_atom.h"

namespace interaction {

class Dihedral_Colvar : public Colvar {
public:
  Dihedral_Colvar() : Colvar("Dihedral") {}

  topology::dihedral_restraint_struct *params;

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
