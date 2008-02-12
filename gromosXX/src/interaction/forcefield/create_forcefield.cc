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

#include <interaction/molecular_virial_interaction.h>

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
  if (create_g96_nonbonded(ff, topo, sim, it, os, quiet))
    return 1;

  // correct the virial (if molecular virial is required)
  if (sim.param().pcouple.virial == math::molecular_virial){
    
    Molecular_Virial_Interaction * mvi = new Molecular_Virial_Interaction;
    ff.push_back(mvi);
  }

  // the special
  DEBUG(8, "creating the special terms");
  if(create_special(ff, topo, sim.param(), os, quiet))
    return 1;

  if (!quiet){
  
    if (sim.param().perturbation.perturbation){
      os << "\t" << std::setw(20) << std::left << "perturbation" 
	 << std::setw(30) << "on" << std::right << "\n"
	 << "\t\tlambda         : " << sim.param().perturbation.lambda << "\n"
	 << "\t\texponent       : " << sim.param().perturbation.lambda_exponent << "\n"
	 << "\t\tdlambda        : " << sim.param().perturbation.dlamt << "\n"
	 << "\t\tscaling        : ";
      
      if (sim.param().perturbation.scaling){
	if (sim.param().perturbation.scaled_only)
	  os << "perturbing only scaled interactions\n";
	else
	  os << "on\n";
      }
      else
	os << "off\n";

      if (topo.perturbed_solute().atoms().size() == 0)
	os << "\t\t" << "using unperturbed nonbonded routines as no atoms are perturbed\n";
      else os << "\t\t" << "with " << topo.perturbed_solute().atoms().size() << " perturbed atoms\n";

    }
    else{
      os << "\t" << std::setw(20) << std::left << "perturbation" 
	 << std::setw(30) << std::left << "off" << std::right << "\n";
    }
    
    os << "END\n";
  }
  
  return 0;

}

