/**
 * @file create_special.cc
 * create the special terms.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <math/periodicity.h>

// special interactions
#include <interaction/interaction_types.h>
#include <interaction/special/position_restraint_interaction.h>

#include <io/instream.h>
#include <io/topology/in_topology.h>

#include "create_special.h"

template<math::virial_enum v>
struct special_interaction_spec
{
  static const math::virial_enum do_virial = v;
};

template<typename t_interaction_spec>
static void _create_special(interaction::Forcefield & ff,
			    topology::Topology const & topo,
			    simulation::Parameter const & param,
			    bool quiet = false)
{
  if (!quiet)
    std::cout << "SPECIAL\n";
  
  if (param.posrest.posrest == 1 || 
      param.posrest.posrest == 2){

    if(!quiet)
      std::cout <<"\tPosition restraints\n";

    interaction::Position_Restraint_Interaction<t_interaction_spec> *pr =
      new interaction::Position_Restraint_Interaction<t_interaction_spec>();

    ff.push_back(pr);
    
    if (param.pcouple.virial == math::atomic_virial)
      io::messages.add("Position restraints with atomic virial ill defined",
		       "create_special", io::message::warning);
    else if (param.pcouple.virial == math::molecular_virial)
      io::messages.add("Position restraint forces not added to molecular virial",
		       "create_special", io::message::warning);
  }
  else if (param.posrest.posrest == 3){
    io::messages.add("Position constraints not implemented",
		     "create_special", io::message::error);
  }
  
}

int interaction::create_special(interaction::Forcefield & ff,
				topology::Topology const & topo,
				simulation::Parameter const & param,
				bool quiet)
{
  switch(param.pcouple.virial){
    case math::no_virial:
      {
	// create an interaction spec suitable for the bonded terms
	_create_special<special_interaction_spec<math::no_virial> >
	  (ff, topo, param, quiet);
	break;
      }
    case math::atomic_virial:
      {
	// create an interaction spec suitable for the bonded terms
	_create_special<special_interaction_spec<math::atomic_virial> >
	  (ff, topo, param, quiet);
	break;
      }
    case math::molecular_virial:
      {
	// create an interaction spec suitable for the bonded terms
	_create_special<special_interaction_spec<math::molecular_virial> >
	  (ff, topo, param, quiet);
	break;
      }
    default:
      {
	throw std::string("Wrong virial type requested");
      }
  }
  return 0;
}

