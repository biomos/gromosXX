/**
 * @file create_forcefield.cc
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
#include <interaction/interaction_types.h>
#include <interaction/forcefield/forcefield.h>

#include <io/instream.h>
#include <io/topology/in_topology.h>


#include "create_forcefield.h"
#include <interaction/bonded/create_bonded.h>
#include <interaction/nonbonded/create_nonbonded.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE forcefield

/**
 * create a Gromos96 (like) forcefield.
 */
void interaction::create_g96_forcefield(interaction::Forcefield & ff,
					topology::Topology const & topo,
					simulation::Parameter const & param,
					io::In_Topology & it)
{
  // the bonded
  DEBUG(8, "creating the bonded terms");
  create_g96_bonded(ff, topo, param, it);

  // the nonbonded
  DEBUG(8, "creating the nonbonded terms");
  create_g96_nonbonded(ff, topo, param, it);
  
  // the special

}

