/**
 * @file create_unperturbed.cc
 * create the non-perturbed nonbonded interaction
 * using a standard pairlist algorithm.
 */

#include <util/stdheader.h>

#ifdef OMP
#include <omp.h>
#endif

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

// nonbonded base
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

// nonbonded pairlist
#include <interaction/nonbonded/pairlist/pairlist.h>


// nonbonded interaction
#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_term.h>
#include <interaction/nonbonded/interaction/nonbonded_innerloop.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_innerloop.h>
#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_pair.h>
#include <interaction/nonbonded/interaction/nonbonded_set.h>
#include <interaction/nonbonded/interaction/nonbonded_interaction.h>

// nonbonded filter
#include <interaction/nonbonded/filter/filter.h>
#include <interaction/nonbonded/filter/exclusion_filter.h>
#include <interaction/nonbonded/filter/chargegroup_grid.h>
#include <interaction/nonbonded/filter/range_filter.h>

// nonbonded pairlist algorithm
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/standard_pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/grid_pairlist_algorithm.h>


#include <interaction/nonbonded/interaction_spec.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

//==================================================
// DECLARATIONS
//==================================================

template<math::boundary_enum t_boundary>
static void _select_virial(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Simulation const & sim,
			   configuration::Configuration const & conf,
			   io::In_Topology & it, bool quiet);

template<math::boundary_enum t_boundary, math::virial_enum t_virial>
static void _select_cutoff(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Simulation const & sim,
			   configuration::Configuration const & conf,
			   io::In_Topology & it, bool quiet);

template<math::boundary_enum t_boundary, 
	 math::virial_enum t_virial,
	 bool t_cutoff>
static void _add_nonbonded(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Simulation const & sim,
			   configuration::Configuration const & conf,
			   io::In_Topology & it, bool quiet);

//==================================================
// DEFINITION
//==================================================

namespace interaction
{
  void create_g96_unperturbed(interaction::Forcefield & ff,
			      topology::Topology const & topo,
			      simulation::Simulation const & sim,
			      configuration::Configuration const & conf,
			      io::In_Topology & it,
			      bool quiet)
  {
    DEBUG(9, "\tcreate g96 nonbonded terms");
    
    if (sim.param().force.nonbonded == 1){
      
      switch(sim.param().boundary.boundary){
	case math::vacuum:
	  {
	    // no pressure calculation in vaccuo
	    _select_cutoff<math::vacuum, math::no_virial>(ff, topo, sim, conf, it, quiet);
	    break;
	  }
	case math::rectangular:
	  {
	  _select_virial<math::rectangular>(ff, topo, sim, conf, it, quiet);
	  break;
	  }
	case math::triclinic:
	  {
	    _select_virial<math::triclinic>(ff, topo, sim, conf, it, quiet);
	    break;
	  }
	default:
	{
	  throw std::string("wrong boundary condition requested");
	}
      }
    }
  }
}

template<math::boundary_enum t_boundary>
static void _select_virial(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Simulation const & sim,
			   configuration::Configuration const & conf,
			   io::In_Topology & it,
			   bool quiet)
{
  DEBUG(9, "\t\tboundary : " << t_boundary);
  
  switch(sim.param().pcouple.virial){
    case math::no_virial:
      {
	_select_cutoff<t_boundary, math::no_virial>(ff, topo, sim, conf, it, quiet);
	break;
      }
    case math::atomic_virial:
      {
	_select_cutoff<t_boundary, math::atomic_virial>(ff, topo, sim, conf, it, quiet);
	break;
      }
    case math::molecular_virial:
      {
	_select_cutoff<t_boundary, math::molecular_virial>(ff, topo, sim, conf, it, quiet);
	break;
      }
    default:
      {
	throw std::string("Wrong virial type requested");
      }
  }
}

template<math::boundary_enum t_boundary, math::virial_enum t_virial>
static void _select_cutoff(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Simulation const & sim,
			   configuration::Configuration const & conf,
			   io::In_Topology & it,
			   bool quiet)
{
  DEBUG(9, "\t\tvirial : " << t_virial);

  if (sim.param().pairlist.atomic_cutoff){
    
    _add_nonbonded<t_boundary, t_virial, interaction::atomic_cutoff_on>
      (ff, topo, sim, conf, it, quiet);
  }
  else
    _add_nonbonded<t_boundary, t_virial, interaction::atomic_cutoff_off>
      (ff, topo, sim, conf, it, quiet);
}

template<math::boundary_enum t_boundary, 
	 math::virial_enum t_virial,
	 bool t_cutoff>
static void _add_nonbonded(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Simulation const & sim,
			   configuration::Configuration const & conf,
			   io::In_Topology & it,
			   bool quiet)
{
  DEBUG(9, "\t\tperturbation : off");
  DEBUG(9, "\t\tscaling      : off");
  DEBUG(9, "\t\tstandard pairlist");
  
  if (!quiet){
    std::cout << "\tnonbonded interaction\n"
	      << "\t\tperturbation      : off\n"
	      << "\t\tscaling           : off\n"
	      << "\t\tvirial            : ";
  
    switch(t_virial){
      case math::no_virial:
	std::cout << "off\n";
	break;
      case math::molecular_virial:
	std::cout << "molecular\n";
	break;
      case math::atomic_virial:
	std::cout << "atomic\n";
      break;
    }

    std::cout << "\t\tcutoff            : ";
    if(t_cutoff)
      std::cout << "atomic\n";
    else
      std::cout << "chargegroup\n";

    if (sim.param().longrange.rf_excluded)
      std::cout << "\t\treaction field contributions from excluded atoms added\n";
    else
      std::cout << "\t\tno reaction field contributions from excluded atoms\n";
    
    std::cout << "\t\tstandard pairlist\n";
    
    std::cout << "\n";
  }
  
  typedef interaction::Interaction_Spec
    < 
    t_boundary,
    t_virial,
    t_cutoff,
    interaction::bekker_off
    >
    interaction_spec_type;

  typedef interaction::Perturbation_Spec
    <
    interaction::perturbation_off,
    interaction::scaling_off
    >
    perturbation_spec_type;
  
  interaction::Standard_Pairlist_Algorithm<interaction_spec_type, perturbation_spec_type>
    * pa 
    = new  interaction::Standard_Pairlist_Algorithm<interaction_spec_type, 
    perturbation_spec_type>;
  

  interaction::Nonbonded_Interaction<interaction_spec_type, perturbation_spec_type> * the_nonbonded 
    = new interaction::Nonbonded_Interaction<interaction_spec_type, perturbation_spec_type>(pa);
  
  it.read_lj_parameter(the_nonbonded->lj_parameter());
  
  the_nonbonded->initialize(topo, conf, sim);
    
  ff.push_back(the_nonbonded);
}
