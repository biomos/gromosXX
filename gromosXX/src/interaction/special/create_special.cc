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

// special interactions
#include <interaction/interaction_types.h>

#include <interaction/special/position_restraint_interaction.h>
#include <interaction/special/jvalue_restraint_interaction.h>

#include <interaction/bonded/dihedral_interaction.h>
#include <interaction/special/pscale.h>

#include <io/instream.h>
#include <io/topology/in_topology.h>

#include "create_special.h"

int interaction::create_special(interaction::Forcefield & ff,
				topology::Topology const & topo,
				simulation::Parameter const & param,
				bool quiet)
{
  if (!quiet)
    std::cout << "SPECIAL\n";
  
  // Position restraints / constraints
  if (param.posrest.posrest == 1 || 
      param.posrest.posrest == 2){

    if(!quiet)
      std::cout <<"\tPosition restraints\n";

    interaction::Position_Restraint_Interaction *pr =
      new interaction::Position_Restraint_Interaction();

    ff.push_back(pr);
    
    if (param.pcouple.virial == math::atomic_virial)
      io::messages.add("Position restraints with atomic virial ill defined",
		       "create_special", io::message::warning);
    else if (param.pcouple.virial == math::molecular_virial)
      io::messages.add("Position restraint forces not added to molecular virial",
		       "create_special", io::message::warning);
  }
  else if (param.posrest.posrest != 0 && param.posrest.posrest != 3){
    io::messages.add("Wrong value for position restraints",
		     "create_special", io::message::error);
  }

  // J-Value restraints
  if (param.jvalue.mode != simulation::restr_off){
    if(!quiet){
      std::cout << "\tJ-Value restraints (";
      switch(param.jvalue.mode){
	case simulation::restr_inst :
	  std::cout << "instantaneous";
	  break;
	case simulation::restr_av :
	  std::cout << "time averaged";
	  break;
	case simulation::restr_biq :
	  std::cout << "biquadratic";
	  break;
	default:
	  std::cout << "unknown mode!";
	  break;
      }
      std::cout << ")\n";
    }

    interaction::Jvalue_Restraint_Interaction *jr =
      new interaction::Jvalue_Restraint_Interaction();
    
    ff.push_back(jr);
  }

  // Periodic Scaling
  // right now this only works if no dihedral angles with j-value restraints on are perturbed...
  if (param.pscale.jrest){

    if(!quiet){
      std::cout << "\tscaling based on J-Value restraints\n";
    }
    
    interaction::Periodic_Scaling * ps = 
      new interaction::Periodic_Scaling(ff, param);

    ff.push_back(ps);
  }

  return 0;
}
