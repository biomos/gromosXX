/**
 * @file create_forcefield.cc
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>
#include <interaction/forcefield/forcefield.h>

#include <io/ifp.h>

#include "create_forcefield.h"
#include <interaction/bonded/create_bonded.h>
#include <interaction/nonbonded/create_nonbonded.h>
#include <interaction/special/create_special.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE forcefield

/**
 * create a Gromos96 (like) forcefield.
 */
int interaction::create_g96_forcefield(interaction::Forcefield & ff,
				       topology::Topology const & topo,
				       simulation::Simulation const & sim,
				       configuration::Configuration const & conf,
				       io::IFP & it,
				       std::ostream & os,
				       bool quiet)
{
  if (!quiet)
    os << "FORCEFIELD\n";
  
  // the bonded
  DEBUG(8, "creating the bonded terms");
  if (create_g96_bonded(ff, topo, sim.param(), it, os, quiet))
    return 1;

  // the nonbonded
  DEBUG(8, "creating the nonbonded terms");
  if (create_g96_nonbonded(ff, topo, sim, conf, it, os, quiet))
    return 1;

  // the special
  DEBUG(8, "creating the special terms");
  if(create_special(ff, topo, sim.param(), os, quiet))
    return 1;

  if (!quiet)
    os << "END\n";

  return 0;

}

