/**
 * @file standard_pairlist_algorithm_atomic.cc
 * standard pairlist algorithm (atomic implementation)
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <math/periodicity.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/nonbonded_innerloop.h>

#include <interaction/nonbonded/interaction/perturbed_nonbonded_term.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_innerloop.h>

#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/standard_pairlist_algorithm.h>

#include <interaction/nonbonded/interaction_spec.h>

#include <util/debug.h>
#include <util/template_split.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist


void interaction::Standard_Pairlist_Algorithm::
update_atomic(topology::Topology & topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim, 
	      interaction::Storage & storage,
	      interaction::Pairlist & pairlist,
	      unsigned int begin, unsigned int end,
	      unsigned int stride)
{
  SPLIT_INNERLOOP(_update_atomic, topo, conf, sim, storage, pairlist, begin, end, stride);
}

void interaction::Standard_Pairlist_Algorithm::
update_perturbed_atomic(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			interaction::Storage & storage,
			interaction::Pairlist & pairlist,
			interaction::Pairlist & perturbed_pairlist,
			unsigned int begin, unsigned int end,
			unsigned int stride)
{
  SPLIT_PERTURBATION(update_pert_atomic,
		     topo, conf, sim, storage,
		     pairlist, perturbed_pairlist,
		     begin, end, stride);
}


template<typename t_interaction_spec>
void interaction::Standard_Pairlist_Algorithm::
_update_atomic(topology::Topology & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation & sim,
	       interaction::Storage & storage,
	       interaction::Pairlist & pairlist,
	       unsigned int begin, unsigned int end,
	       unsigned int stride)
{
  DEBUG(7, "standard pairlist update (atomic cutoff)");
  const double update_start = util::now();
  
  // create the innerloop
  Nonbonded_Innerloop innerloop(*m_param);
  innerloop.init(sim);

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  math::VArray const & pos = conf.current().pos;
  math::Vec v;

  // empty the pairlist
  for(unsigned int i=0; i<topo.num_atoms(); ++i)
    pairlist[i].clear();

  DEBUG(7, "pairlist empty");

  const int num_solute = topo.num_solute_atoms();
  const int num_atoms = topo.num_atoms();
  
  int a1 = begin;

  for( ; a1 < num_solute; a1 += stride) {

    DEBUG(9, "solute (" << a1 << ") - solute");
    
    for(int a2 = a1 + 1; a2 < num_solute; ++a2){

      assert(a1 != a2);
      
      // check exclusions and range
      periodicity.nearest_image(pos(a1), pos(a2), v);
      
      // the distance
      const double d = math::abs2(v);

      DEBUG(10, "\t" << a1 << " - " << a2);
    
      if (d > m_cutoff_long_2){        // OUTSIDE
	DEBUG(11, "\t\toutside");
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	DEBUG(11, "\t\tlongrange");
	// the interactions
	innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);

	continue;
      } // longrange
      
      // shortrange - check exclusions
      if (excluded_solute_pair(topo, a1, a2)){
	continue;
      }

      DEBUG(11, "\t\tshortrange");
      pairlist[a1].push_back(a2);
      
    } // solute - solute

    DEBUG(9, "solute (" << a1 << ") - solvent");

    for(int a2 = num_solute; a2 < num_atoms; ++a2){
    
      assert(a1 != a2);
      
      // check range, no exclusions
      periodicity.nearest_image(pos(a1), pos(a2), v);
      
      // the distance
      const double d = math::abs2(v);
    
      DEBUG(10, "\t" << a1 << " - " << a2);

      if (d > m_cutoff_long_2){        // OUTSIDE
	DEBUG(11, "\t\toutside");
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	DEBUG(11, "\t\tlongrange");
	// the interactions
	innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);

	continue;
      } // longrange
      
      DEBUG(11, "\t\tshortrange");
      pairlist[a1].push_back(a2);

    } // solute - solvent
    
  }

  // int a1 = num_solute;

  // multiple solvents
  DEBUG(9, "solvent - solvent");

  for(unsigned int s = 0; s < topo.num_solvents(); ++s){
    DEBUG(11, "solvent " << s);
    int end = a1 + topo.num_solvent_molecules(s);
    DEBUG(11, "\tends at atom " << end);
    
    const int num_solv_at = topo.num_solvent_atoms(s);
    int a2_start = a1 + num_solv_at;
      
    DEBUG(11, "\twith " << num_solv_at << " atoms");
    
    for( ; a1 < end; a1 += stride){
      
      if (a1 == a2_start) a2_start += num_solv_at;
      
      for(int a2 = a2_start; a2_start < num_atoms; ++a2){
	
	assert(a1 != a2);
	
	// check range, no exclusions
	periodicity.nearest_image(pos(a1), pos(a2), v);
	
	// the distance
	const double d = math::abs2(v);
	
	DEBUG(10, "\t" << a1 << " - " << a2);
	
	if (d > m_cutoff_long_2){        // OUTSIDE
	  DEBUG(11, "\t\toutside");
	  continue;
	}
	
	if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	  DEBUG(11, "\t\tlongrange");	  
	  // the interactions
	  innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);
	  
	  continue;
	} // longrange
	
	  DEBUG(11, "\t\tshortrange");	
	  pairlist[a1].push_back(a2);
	  
      } // solvent - solvent
      
    } // a1 of solvent s
    
  } // multiple solvents
  
  this->m_timing += util::now() - update_start;
  
  DEBUG(7, "pairlist done");

}


template<typename t_perturbation_details>
void interaction::Standard_Pairlist_Algorithm::
update_pert_atomic(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim, 
		   interaction::Storage & storage,
		   interaction::Pairlist & pairlist,
		   interaction::Pairlist & perturbed_pairlist,
		   unsigned int begin, unsigned int end,
		   unsigned int stride)
{
  SPLIT_PERT_INNERLOOP(_update_pert_atomic,
		       topo, conf, sim, storage,
		       pairlist, perturbed_pairlist, 
		       begin, end, stride);
}

template<typename t_interaction_spec, typename t_perturbation_details>
void interaction::Standard_Pairlist_Algorithm::
_update_pert_atomic(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    interaction::Storage & storage,
		    interaction::Pairlist & pairlist,
		    interaction::Pairlist & perturbed_pairlist,
		    unsigned int begin, unsigned int end,
		    unsigned int stride)
{
  DEBUG(7, "standard pairlist update (atomic cutoff)");
  const double update_start = util::now();
  
  // create the innerloops
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  
  Nonbonded_Innerloop innerloop(*m_param);
  innerloop.init(sim);
  
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details>
    perturbed_innerloop(*m_param);
  perturbed_innerloop.init(sim);
  perturbed_innerloop.set_lambda(topo.lambda(), topo.lambda_exp());

  // empty the pairlist
  assert(pairlist.size() == topo.num_atoms());
  assert(perturbed_pairlist.size() == topo.num_atoms());
  
  for(unsigned int i=0; i<topo.num_atoms(); ++i){
    pairlist[i].clear();
    perturbed_pairlist[i].clear();
  }

  DEBUG(7, "pairlist cleared");

  const int num_solute = topo.num_solute_atoms();
  const int num_atoms = topo.num_atoms();

  math::VArray const & pos = conf.current().pos;
  math::Vec v;

  int a1 = begin;

  for( ; a1 < num_solute; a1 += stride) {
    DEBUG(9, "solute (" << a1 << ") - solute");

    for(int a2 = a1 + 1; a2 < num_solute; ++a2){

      assert(a1 != a2);
      
      // check exclusions and range
      periodicity.nearest_image(pos(a1), pos(a2), v);
      
      // the distance
      const double d = math::abs2(v);

      DEBUG(10, "\t" << a1 << " - " << a2);
    
      if (d > m_cutoff_long_2){        // OUTSIDE
	DEBUG(11, "\t\toutside");
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	DEBUG(11, "\t\tlongrange");
	
	// the interactions
	if (topo.is_perturbed(a1)){
	  DEBUG(11, "\t\t" << a1 << " perturbed");

	  if (t_perturbation_details::do_scaling &&
	      sim.param().perturbation.scaled_only){

	    // ok, only perturbation if it is a scaled pair...
	    std::pair<int, int> 
	      energy_group_pair(topo.atom_energy_group(a1),
				topo.atom_energy_group(a2));
	    
	    if (topo.energy_group_scaling().count(energy_group_pair))
	      perturbed_innerloop.
		perturbed_lj_crf_innerloop(topo, conf, a1, a2,
					   storage, periodicity);
	    else
	      innerloop.
		lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2,
						     storage, periodicity);
	  } // perturb scaled interactions only
	  else{
	    perturbed_innerloop.
	      perturbed_lj_crf_innerloop(topo, conf, a1, a2,
					 storage, periodicity);
	  }
	}
	else if (topo.is_perturbed(a2)){
	  DEBUG(11, "\t\t" << a2 << " perturbed");

	  if (t_perturbation_details::do_scaling &&
	      sim.param().perturbation.scaled_only){

	    // ok, only perturbation if it is a scaled pair...
	    std::pair<int, int> 
	      energy_group_pair(topo.atom_energy_group(a1),
				topo.atom_energy_group(a2));
	    
	    if (topo.energy_group_scaling().count(energy_group_pair))
	      perturbed_innerloop.
		perturbed_lj_crf_innerloop(topo, conf, a2, a1,
					   storage, periodicity);
	    else
	      innerloop.
		lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2,
						     storage, periodicity);
	  } // perturb scaled interactions only
	  else{
	    perturbed_innerloop.
	      perturbed_lj_crf_innerloop(topo, conf, a2, a1,
					 storage, periodicity);
	  }
	}
	else{
	  DEBUG(11, "\t\tnot perturbed");
	  innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);
	}
	
	continue;
      } // longrange
      
      // shortrange - check exclusions
      if (excluded_solute_pair(topo, a1, a2)) continue;

      DEBUG(11, "\t\tshortrange");

      if (topo.is_perturbed(a1)){
	DEBUG(11, "\t\t" << a1 << " perturbed");

	if (t_perturbation_details::do_scaling &&
	    sim.param().perturbation.scaled_only){

	  // ok, only perturbation if it is a scaled pair...
	  std::pair<int, int> 
	    energy_group_pair(topo.atom_energy_group(a1),
			      topo.atom_energy_group(a2));
	  
	  if (topo.energy_group_scaling().count(energy_group_pair))
	    perturbed_pairlist[a1].push_back(a2);
	  else
	    pairlist[a1].push_back(a2);
	} // scaling
	else{
	  perturbed_pairlist[a1].push_back(a2);
	}
      }
      else if (topo.is_perturbed(a2)){
	DEBUG(11, "\t\t" << a2 << " perturbed");

	if (t_perturbation_details::do_scaling &&
	    sim.param().perturbation.scaled_only){

	  // ok, only perturbation if it is a scaled pair...
	  std::pair<int, int> 
	    energy_group_pair(topo.atom_energy_group(a1),
			      topo.atom_energy_group(a2));
	  
	  if (topo.energy_group_scaling().count(energy_group_pair))
	    perturbed_pairlist[a2].push_back(a1);
	  else
	    pairlist[a1].push_back(a2);
	} // scaling
	else{
	  perturbed_pairlist[a2].push_back(a1);
	}
	
      }
      else{
	DEBUG(11, "\t\tnot perturbed");
	pairlist[a1].push_back(a2);
      }
      
    } // solute - solute

    DEBUG(9, "solute (" << a1 << ") - solvent");

    for(int a2 = num_solute; a2 < num_atoms; ++a2){
    
      assert(a1 != a2);
      
      // check range, no exclusions
      periodicity.nearest_image(pos(a1), pos(a2), v);
      
      // the distance
      const double d = math::abs2(v);
    
      DEBUG(10, "\t" << a1 << " - " << a2);

      if (d > m_cutoff_long_2){        // OUTSIDE
	DEBUG(11, "\t\toutside");
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	DEBUG(11, "\t\tlongrange");
	
	// the interactions
	if (topo.is_perturbed(a1)){

	  DEBUG(11, "\t\t" << a1 << " perturbed");
	  if (t_perturbation_details::do_scaling &&
	      sim.param().perturbation.scaled_only){

	    // ok, only perturbation if it is a scaled pair...
	    std::pair<int, int> 
	      energy_group_pair(topo.atom_energy_group(a1),
				topo.atom_energy_group(a2));
	    
	    if (topo.energy_group_scaling().count(energy_group_pair))
	      perturbed_innerloop.
		perturbed_lj_crf_innerloop(topo, conf, a1, a2, storage, periodicity);
	    else
	      innerloop.
		lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);

	  } // scaling
	  else{
	    perturbed_innerloop.
	      perturbed_lj_crf_innerloop(topo, conf, a1, a2, storage, periodicity);
	  }
	}
	else{
	  DEBUG(11, "\t\tnot perturbed");
	  innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);
	}
	
	continue;
      } // longrange
      
      DEBUG(11, "\t\tshortrange");

      if (topo.is_perturbed(a1)){
	DEBUG(11, "\t\t" << a1 << " perturbed");

	if (t_perturbation_details::do_scaling &&
	    sim.param().perturbation.scaled_only){

	  // ok, only perturbation if it is a scaled pair...
	  std::pair<int, int> 
	    energy_group_pair(topo.atom_energy_group(a1),
			      topo.atom_energy_group(a2));
	  
	  if (topo.energy_group_scaling().count(energy_group_pair))
	    perturbed_pairlist[a1].push_back(a2);
	  else
	    pairlist[a1].push_back(a2);

	} // scaling
	else{
	  perturbed_pairlist[a1].push_back(a2);
	}
      }
      else{
	DEBUG(11, "\t\tnot perturbed");
	pairlist[a1].push_back(a2);
      }

    } // solute - solvent
    
  }

  int solv_start = num_solute;

  // multiple solvents
  DEBUG(9, "solvent - solvent");
  DEBUG(10, "\tnum_atoms = " << num_atoms);

  for(unsigned int s = 0; s < topo.num_solvents(); ++s){
    DEBUG(11, "solvent " << s);

    int end = solv_start + topo.num_solvent_atoms(s);
    DEBUG(11, "\tends at atom " << end);

    const int num_solv_at = topo.num_solvent_atoms(s) / topo.num_solvent_molecules(s);
    int a2_start = solv_start + num_solv_at;
    DEBUG(11, "\twith " << num_solv_at << " atoms");
    DEBUG(11, "\ta1 starts with " << a1 << "\ta2 starts with " << a2_start);

    for( ; a1 < end; a1+=stride){
      
      while (a1 >= a2_start) a2_start += num_solv_at;
      
      for(int a2 = a2_start; a2 < num_atoms; ++a2){
	
	assert(a1 != a2);
	
	// check range, no exclusions
	periodicity.nearest_image(pos(a1), pos(a2), v);
	
	// the distance
	const double d = math::abs2(v);

	DEBUG(10, "\t" << a1 << " - " << a2);
	
	if (d > m_cutoff_long_2){        // OUTSIDE
	  DEBUG(11, "\t\toutside");
	  continue;
	}
  
	if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	  DEBUG(11, "\t\tlongrange");	  

	  // the interactions
	  innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);
	  
	  continue;
	} // longrange

	DEBUG(11, "\t\tshortrange");	
	pairlist[a1].push_back(a2);
	
      } // solvent - solvent

    } // a1 of solvent s

    // start of next solvent
    solv_start += topo.num_solvent_atoms(s);    

  } // multiple solvents
  

  this->m_timing += util::now() - update_start;
  
  DEBUG(7, "pairlist done");

}
