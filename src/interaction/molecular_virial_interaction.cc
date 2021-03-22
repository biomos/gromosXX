/**
 * @file molecular_virial_interaction.cc
 * recover molecular virial from atomic virial
 */

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../interaction/interaction.h"

#include "molecular_virial_interaction.h"

#include "../util/prepare_virial.h"

#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

int interaction::Molecular_Virial_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();
  util::atomic_to_molecular_virial(topo, conf, sim);
  m_timer.reset();
  return 0;
  
}
