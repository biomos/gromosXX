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

// nonbonded base
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_base.h>

// nonbonded filter
#include <interaction/nonbonded/filter/filter.h>
#include <interaction/nonbonded/filter/exclusion_filter.h>
#include <interaction/nonbonded/filter/perturbation_filter.h>
#include <interaction/nonbonded/filter/chargegroup_grid.h>
#include <interaction/nonbonded/filter/range_filter.h>

// nonbonded pairlist algorithm
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/standard_pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/grid_pairlist_algorithm.h>

// nonbonded pairlist
#include <interaction/nonbonded/pairlist/pairlist.h>


// nonbonded interaction
#include <interaction/nonbonded/interaction/nonbonded_innerloop.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_innerloop.h>
#include <interaction/nonbonded/interaction/nonbonded_interaction.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_interaction.h>

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
			   simulation::Parameter const & param,
			   io::In_Topology & it);

template<math::boundary_enum t_boundary, math::virial_enum t_virial>
static void _select_cutoff(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Parameter const & param,
			   io::In_Topology & it);

template<math::boundary_enum t_boundary, 
	 math::virial_enum t_virial,
	 bool t_cutoff>
static void _add_grid_nonbonded(interaction::Forcefield & ff,
				topology::Topology const & topo,
				simulation::Parameter const & param,
				io::In_Topology & it);

//==================================================
// DEFINITION
//==================================================

namespace interaction
{
  
  void create_g96_unperturbed_grid(interaction::Forcefield & ff,
				   topology::Topology const & topo,
				   simulation::Parameter const & param,
				   io::In_Topology & it)
  {
    DEBUG(9, "\tcreate g96 nonbonded terms");
    
    if (param.force.nonbonded == 1){
      
      switch(param.boundary.boundary){
	case math::vacuum:
	  {
	  // no pressure calculation in vaccuo
	    _select_cutoff<math::vacuum, math::no_virial>(ff, topo, param, it);
	    break;
	  }
	case math::rectangular:
	  {
	    _select_virial<math::rectangular>(ff, topo, param, it);
	    break;
	  }
	case math::triclinic:
	  {
	    _select_virial<math::triclinic>(ff, topo, param, it);
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
			   simulation::Parameter const & param,
			   io::In_Topology & it)
{
  DEBUG(9, "\t\tboundary : " << t_boundary);
  
  switch(param.pcouple.virial){
    case math::no_virial:
      {
	_select_cutoff<t_boundary, math::no_virial>(ff, topo, param, it);
	break;
      }
    case math::atomic_virial:
      {
	_select_cutoff<t_boundary, math::atomic_virial>(ff, topo, param, it);
	break;
      }
    case math::molecular_virial:
      {
	_select_cutoff<t_boundary, math::molecular_virial>(ff, topo, param, it);
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
			   simulation::Parameter const & param,
			   io::In_Topology & it)
{
  DEBUG(9, "\t\tvirial : " << t_virial);

  if (param.pairlist.atomic_cutoff){
    
    _add_grid_nonbonded<t_boundary, t_virial, true>(ff, topo, param, it);
  }
  else
    _add_grid_nonbonded<t_boundary, t_virial, false>(ff, topo, param, it);
}

template<math::boundary_enum t_boundary, 
	 math::virial_enum t_virial,
	 bool t_cutoff>
static void _add_grid_nonbonded(interaction::Forcefield & ff,
				topology::Topology const & topo,
				simulation::Parameter const & param,
				io::In_Topology & it)
{
  DEBUG(9, "\t\tperturbation : off");
  DEBUG(9, "\t\tscaling      : off");
  DEBUG(9, "\t\tgrid based pairlist");

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

  if (param.longrange.rf_excluded)
    std::cout << "\t\treaction field contributions from excluded atoms added\n";
  else
    std::cout << "\t\tno reaction field contributions from excluded atoms\n";

  std::cout << "\t\tgrid based pairlist\n";

  std::cout << "\n";

  typedef interaction::Nonbonded_Interaction
    < 
    interaction::Grid_Interaction_Spec
    < 
    t_boundary,
    false,
    t_virial,
    t_cutoff,
    false
    >
    >
    nonbonded_type;
    
  nonbonded_type * the_nonbonded = new nonbonded_type;
  
  it.read_lj_parameter(the_nonbonded->lj_parameter());
  
  the_nonbonded->energies.resize(param.force.energy_group.size(),
				 param.multibath.multibath.size());

  the_nonbonded->perturbed_energy_derivatives.
    resize(param.force.energy_group.size(),
	   param.multibath.multibath.size());
  
  ff.push_back(the_nonbonded);
}

