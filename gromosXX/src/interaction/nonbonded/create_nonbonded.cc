/**
 * @file create_nonbonded.cc
 * create the nonbonded interaction.
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
			      simulation::Simulation const & sim,
			      configuration::Configuration const & conf,
			      io::In_Topology & it);
  
  void create_g96_unperturbed_grid(interaction::Forcefield & ff,
				   topology::Topology const & topo,
				   simulation::Simulation const & sim,
				   configuration::Configuration const & conf,
				   io::In_Topology & it);
  
  void create_g96_perturbed(interaction::Forcefield & ff,
			    topology::Topology const & topo,
			    simulation::Simulation const & sim,
			    configuration::Configuration const & conf,
			    io::In_Topology & it);
  
  void create_g96_perturbed_grid(interaction::Forcefield & ff,
				 topology::Topology const & topo,
				 simulation::Simulation const & sim,
				 configuration::Configuration const & conf,
				 io::In_Topology & it);
}

//==================================================
// DEFINITION
//==================================================

void interaction::create_g96_nonbonded(interaction::Forcefield & ff,
				       topology::Topology const & topo,
				       simulation::Simulation const & sim,
				       configuration::Configuration const & conf,
				       io::In_Topology & it)
{
  DEBUG(9, "\tcreate g96 nonbonded terms");

  if (sim.param().force.nonbonded == 1){

    if (sim.param().perturbation.perturbation){
      
      if (sim.param().pairlist.grid){
	create_g96_perturbed_grid(ff, topo, sim,conf, it);
	std::cout << "\t\t\tgrid size          : " << sim.param().pairlist.grid_size << "\n";
      }
      else
	create_g96_perturbed(ff, topo, sim, conf, it);
    }
    else{

      if (sim.param().pairlist.grid){
	create_g96_unperturbed_grid(ff, topo, sim, conf, it);
	std::cout << "\t\tgrid size :     " << sim.param().pairlist.grid_size << "\n";
      }
      else
	create_g96_unperturbed(ff, topo, sim, conf, it);

    }

    std::cout << "\t\t\tinner cutoff           : " << sim.param().pairlist.cutoff_short << "\n"
	      << "\t\t\touter cutoff           : " << sim.param().pairlist.cutoff_long << "\n"
	      << "\t\t\tepsilon                : " << sim.param().longrange.epsilon << "\n"
	      << "\t\t\treactionfield epsilon  : " << sim.param().longrange.rf_epsilon << "\n"
	      << "\t\t\tkappa                  : " << sim.param().longrange.rf_kappa << "\n"
	      << "\t\t\treactionfield cutoff   : " << sim.param().longrange.rf_cutoff << "\n"
	      << "\t\t\tpairlist creation every  " << sim.param().pairlist.skip_step << " steps\n"
	      << "\n";
    
  }
}
