/**
 * @file create_nonbonded.h
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
#include <interaction/forcefield/forcefield.h>

#include <interaction/interaction_types.h>
#include <math/periodicity.h>

#include <io/instream.h>
#include <io/topology/in_topology.h>

#include "create_nonbonded.h"

// general
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

//==================================================
// DECLARATIONS
//==================================================

namespace interaction
{
  void create_g96_unperturbed(interaction::Forcefield & ff,
			      topology::Topology const & topo,
			      simulation::Parameter const & param,
			      io::In_Topology & it);
  
  void create_g96_unperturbed_grid(interaction::Forcefield & ff,
				   topology::Topology const & topo,
				   simulation::Parameter const & param,
				   io::In_Topology & it);
  
  void create_g96_perturbed(interaction::Forcefield & ff,
			    topology::Topology const & topo,
			    simulation::Parameter const & param,
			    io::In_Topology & it);
  
  void create_g96_perturbed_grid(interaction::Forcefield & ff,
				 topology::Topology const & topo,
				 simulation::Parameter const & param,
				 io::In_Topology & it);
}

//==================================================
// DEFINITION
//==================================================

void interaction::create_g96_nonbonded(interaction::Forcefield & ff,
				       topology::Topology const & topo,
				       simulation::Parameter const & param,
				       io::In_Topology & it)
{
  DEBUG(9, "\tcreate g96 nonbonded terms");

  if (param.force.nonbonded == 1){

    if (param.perturbation.perturbation){
      
      if (param.pairlist.grid){
	create_g96_perturbed_grid(ff, topo, param, it);
	std::cout << "\t\t\tgrid size          : " << param.pairlist.grid_size << "\n";
      }
      else
	create_g96_perturbed(ff, topo, param, it);
    }
    else{

      if (param.pairlist.grid){
	create_g96_unperturbed_grid(ff, topo, param, it);
	std::cout << "\t\tgrid size :     " << param.pairlist.grid_size << "\n";
      }
      else
	create_g96_unperturbed(ff, topo, param, it);

    }

    std::cout << "\t\t\tinner cutoff           : " << param.pairlist.cutoff_short << "\n"
	      << "\t\t\touter cutoff           : " << param.pairlist.cutoff_long << "\n"
	      << "\t\t\tepsilon                : " << param.longrange.epsilon << "\n"
	      << "\t\t\treactionfield epsilon  : " << param.longrange.rf_epsilon << "\n"
	      << "\t\t\tkappa                  : " << param.longrange.rf_kappa << "\n"
	      << "\t\t\treactionfield cutoff   : " << param.longrange.rf_cutoff << "\n"
	      << "\t\t\tpairlist creation every  " << param.pairlist.skip_step << " steps\n"
	      << "\n";
    
  }
}
