/**
 * @file external_interaction.cc
 * methods of External_Interaction
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../util/virtual_grain.h"

#include "../../interaction/special/external_interaction.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * add forces from virtual atoms (grains)
 */
int interaction::External_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  // virtual graining external force
  if (cg_topo != NULL &&
      cg_conf != NULL){
    util::update_virtual_force(*cg_topo, *cg_conf, topo, conf, sim);
  }
  
  return 0;
}
