/**
 * @file forcefield.tcc
 * contains the inline functions for
 * forcefield.
 */

#include <util/stdheader.h>

#include <topology/core/core.h>
#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm.h>
#include <interaction/interaction.h>

#include "forcefield.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE forcefield

interaction::Forcefield::~Forcefield()
{
  for(iterator it = begin(), to = end();
      it != to;
      ++it){
    delete *it;
  }
}

int interaction::Forcefield
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim)
{
  DEBUG(5, "forcefield: calculate interaction");

  conf.current().force = 0.0;

  DEBUG(15, "zero energies");
  conf.current().energies.zero();
  DEBUG(15, "zero lambda energies");
  conf.current().perturbed_energy_derivatives.zero();
  conf.current().virial_tensor = 0.0;

  /*
  // prepare for the virial
  if (t_interaction_spec::do_virial == interaction::molecular_virial){
    if(sim.pressure_calculation())
      sim.calculate_mol_com();
  }
  if (t_interaction_spec::do_virial == interaction::atomic_virial){
    if (sim.pressure_calculation())
      sim.calculate_atom_ekin();
  }
  */

  for(iterator it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(7, "interaction: " << (*it)->name);
    (*it)->calculate_interactions(topo, conf, sim);
  }

  return 0;
}
