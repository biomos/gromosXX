/**
 * @file create_forcefield.cc
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
#include <interaction/interaction_types.h>
#include <interaction/forcefield/forcefield.h>

#include <io/instream.h>
#include <io/topology/in_topology.h>


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
				       io::In_Topology & it, bool quiet)
{
  if (!quiet)
    std::cout << "FORCEFIELD\n";
  
  // the bonded
  DEBUG(8, "creating the bonded terms");
  if (create_g96_bonded(ff, topo, sim.param(), it, quiet))
    return 1;

  // the nonbonded
  DEBUG(8, "creating the nonbonded terms");
  if (create_g96_nonbonded(ff, topo, sim, conf, it, quiet))
    return 1;

  // the special
  DEBUG(8, "creating the special terms");
  if(create_special(ff, topo, sim.param(), quiet))
    return 1;

  if (!quiet)
    std::cout << "END\n";

  return 0;

}

