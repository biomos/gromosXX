/**
 * @file standard_pairlist_algorithm_atomic.cc
 * standard pairlist algorithm (atomic implementation)
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction_types.h"
#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"

#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../../../interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"

#include "../../../util/debug.h"
#include "../../../util/template_split.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist


void interaction::Standard_Pairlist_Algorithm::
update_atomic(topology::Topology & topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim, 
	      interaction::PairlistContainer & pairlist,
	      unsigned int begin, unsigned int end,
	      unsigned int stride)
{
  SPLIT_BOUNDARY(_update_atomic, topo, conf, sim, pairlist, begin, end, stride);
}

void interaction::Standard_Pairlist_Algorithm::
update_perturbed_atomic(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			interaction::PairlistContainer & pairlist,
			interaction::PairlistContainer & perturbed_pairlist,
			unsigned int begin, unsigned int end,
			unsigned int stride)
{
  SPLIT_BOUNDARY(_update_pert_atomic,
		       topo, conf, sim,
		       pairlist, perturbed_pairlist, 
		       begin, end, stride);
}


template<math::boundary_enum b>
void interaction::Standard_Pairlist_Algorithm::
_update_atomic(topology::Topology & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation & sim,
	       interaction::PairlistContainer & pairlist,
	       unsigned int begin, unsigned int end,
	       unsigned int stride)
{
  DEBUG(7, "standard pairlist update (atomic cutoff)");
  timer().start_subtimer("pairlist");

  math::Periodicity<b> periodicity(conf.current().box);
  math::VArray const & pos = conf.current().pos;
  math::Vec v;

  // empty the pairlist
  pairlist.clear();

  DEBUG(7, "pairlist empty");

  const int num_solute = topo.num_solute_atoms();
  const int num_atoms = topo.num_atoms();
  
  const simulation::qmmm_enum qmmm = sim.param().qmmm.qmmm;
  
  int a1 = begin;
  DEBUG(9, "\tbegin (1st) = " << begin);  
  for( ; a1 < num_solute; a1 += stride) {

    DEBUG(9, "solute (" << a1 << ") - solute");
    // If a1 is QM
    if (Pairlist_Algorithm::qm_excluded(topo, qmmm, a1)) {
      DEBUG(9, "Skipping all: " << a1 );
      continue;
    }
    for(int a2 = a1 + 1; a2 < num_solute; ++a2){

      assert(a1 != a2);
      // If a2 is QM
      if (Pairlist_Algorithm::qm_excluded(topo, qmmm, a1, a2)) {
			  DEBUG(9, "Skipping pair: " << a1 << "-" << a2);
        continue;
      }
      // check exclusions and range
      periodicity.nearest_image(pos(a1), pos(a2), v);

      // the distance
      const double d = math::abs2(v);

      DEBUG(10, "\t" << a1 << " - " << a2);
    
      if (d > m_cutoff_long_2){        // OUTSIDE
	DEBUG(11, "\t\toutside");
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE
	DEBUG(11, "\t\tlongrange");
        pairlist.solute_long[a1].push_back(a2);

	continue;
      } // longrange
      
      // shortrange - check exclusions
      if (excluded_solute_pair(topo, a1, a2)){
	continue;
      }

      DEBUG(11, "\t\tshortrange");
      pairlist.solute_short[a1].push_back(a2);
      
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
  
      if (d > m_cutoff_short_2){       // LONGRANGE
	DEBUG(11, "\t\tlongrange");
	pairlist.solute_long[a1].push_back(a2);

	continue;
      } // longrange
      
      DEBUG(11, "\t\tshortrange");
      pairlist.solute_short[a1].push_back(a2);

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
    
    if (topo.num_solvent_molecules(s) == 0)
      continue;

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

	DEBUG(11, "\t" << a1 << " - " << a2);
	
	if (d > m_cutoff_long_2){        // OUTSIDE
	  DEBUG(11, "\t\toutside");
	  continue;
	}
  
	if (d > m_cutoff_short_2){       // LONGRANGE
	  DEBUG(11, "\t\tlongrange");	  

          pairlist.solvent_long[a1].push_back(a2);
	  continue;
	} // longrange

	DEBUG(11, "\t\tshortrange");	
	pairlist.solvent_short[a1].push_back(a2);
	
      } // solvent - solvent

    } // a1 of solvent s

    // start of next solvent
    solv_start += topo.num_solvent_atoms(s);    

  } // multiple solvents
  timer().stop_subtimer("pairlist");
  
  DEBUG(7, "pairlist done");

}

template<math::boundary_enum b>
void interaction::Standard_Pairlist_Algorithm::
_update_pert_atomic(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    interaction::PairlistContainer & pairlist,
		    interaction::PairlistContainer & perturbed_pairlist,
		    unsigned int begin, unsigned int end,
		    unsigned int stride)
{
  DEBUG(7, "perturbed standard pairlist update (atomic cutoff)");
  timer().start_subtimer("perturbed pairlist");
  
  math::Periodicity<b> periodicity(conf.current().box);

  // check whether we do scaling && scaling only
  bool scaled_only = (sim.param().perturbation.scaling && sim.param().perturbation.scaled_only);
  
  // empty the pairlist
  assert(pairlist.size() == topo.num_atoms());
  assert(perturbed_pairlist.size() == topo.num_atoms());
  
  pairlist.clear();
  perturbed_pairlist.clear();

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
  
      if (d > m_cutoff_short_2){       // LONGRANGE
	DEBUG(11, "\t\tlongrange");
        
        if (insert_pair(topo, pairlist.solute_long, perturbed_pairlist.solute_long,
                a1, a2, scaled_only))
        {}
        else if (insert_pair(topo, pairlist.solute_long, perturbed_pairlist.solute_long,
                a2, a1, scaled_only))
        {}
        else
          pairlist.solute_long[a1].push_back(a2);
	
	continue;
      } // longrange
      
      // shortrange - check exclusions
      if (excluded_solute_pair(topo, a1, a2)) continue;

      DEBUG(11, "\t\tshortrange");
      
      if (insert_pair(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
              a1, a2, scaled_only))
      {}
      else if (insert_pair(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
              a2, a1, scaled_only))
      {}
      else
        pairlist.solute_short[a1].push_back(a2);

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
  
      if (d > m_cutoff_short_2){       // LONGRANGE
	DEBUG(11, "\t\tlongrange");
	
        if (insert_pair(topo, pairlist.solute_long, perturbed_pairlist.solute_long,
                a1, a2, scaled_only))
        {}
        else
          pairlist.solute_long[a1].push_back(a2);
        
        continue;
      } // longrange
      
      DEBUG(11, "\t\tshortrange");
      
      if (insert_pair(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
              a1, a2, scaled_only))
        ;
      else
        pairlist.solute_short[a1].push_back(a2);

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

    if (topo.num_solvent_molecules(s) == 0)
      continue;

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
  
	if (d > m_cutoff_short_2){       // LONGRANGE
	  DEBUG(11, "\t\tlongrange");	  

          pairlist.solvent_long[a1].push_back(a2);
	  continue;
	} // longrange

	DEBUG(11, "\t\tshortrange");	
	pairlist.solvent_short[a1].push_back(a2);
	
      } // solvent - solvent

    } // a1 of solvent s

    // start of next solvent
    solv_start += topo.num_solvent_atoms(s);    

  } // multiple solvents
  
  timer().stop_subtimer("perturbed pairlist");
  
  DEBUG(7, "pairlist done");

}

// this function has to be doubled because of inlining
inline bool interaction::Standard_Pairlist_Algorithm::insert_pair
(
 topology::Topology & topo,
 interaction::Pairlist & pairlist,
 interaction::Pairlist & perturbed_pairlist,
 int a1, int a2,
 bool scaled_only
 )
{
  if (topo.is_perturbed(a1) || topo.is_eds_perturbed(a1)){
    if (scaled_only){
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
    return true;
  }

  return false;
}

