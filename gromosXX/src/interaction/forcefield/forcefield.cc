/**
 * @file forcefield.cc
 * contains the inline functions for
 * forcefield.
 */

#include <util/stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm.h>
#include <interaction/interaction.h>

#include <util/prepare_virial.h>

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

  const double start = util::now();

  conf.current().force = 0.0;

  DEBUG(15, "zero energies");
  conf.current().energies.zero();
  DEBUG(15, "zero lambda energies");
  conf.current().perturbed_energy_derivatives.zero();
  conf.current().virial_tensor = 0.0;

  // prepare for the virial
  util::prepare_virial(topo, conf, sim);
  
  for(iterator it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(7, "interaction: " << (*it)->name);
    (*it)->calculate_interactions(topo, conf, sim);
  }

  m_timing += util::now() - start;

  return 0;
}


void interaction::Forcefield
::print_timing(std::ostream & os)
{
  os << "    " 
     << std::setw(40) << std::left << name
     << std::setw(20) << m_timing << "\n";

  for(iterator it = begin(), to = end();
      it != to;
      ++it){

    (*it)->print_timing(os);

  }
}

