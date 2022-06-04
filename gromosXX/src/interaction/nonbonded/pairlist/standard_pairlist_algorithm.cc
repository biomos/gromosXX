/**
 * @file standard_pairlist_algorithm.cc
 * standard pairlist algorithm
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

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

interaction::Standard_Pairlist_Algorithm::
Standard_Pairlist_Algorithm()
: interaction::Pairlist_Algorithm(),
  m_solvent_solvent_timing(0.0) {}
        
 /**
   * put the chargegroups into the box
   */
template<math::boundary_enum b>
void _prepare_cog(configuration::Configuration & conf,
topology::Topology & topo)
{
  DEBUG(10, "putting chargegroups into box");
  math::Periodicity<b> periodicity(conf.current().box);
  periodicity.put_chargegroups_into_box(conf, topo);
}

/**
 * calculate center of geometries
 */
int interaction::Standard_Pairlist_Algorithm::
prepare(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "standard pairlist algorithm : prepare");
  
  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);
  
  if (!sim.param().pairlist.atomic_cutoff){

    // first put the chargegroups into the box
    SPLIT_BOUNDARY(_prepare_cog, conf, topo);

    // calculate cg cog's
    DEBUG(10, "calculating cg cog (" << topo.num_solute_chargegroups() << ")");
    m_cg_cog.resize(topo.num_solute_chargegroups());
    math::VArray const &pos = conf.current().pos;
    DEBUG(10, "pos.size() = " << pos.size());

    // calculate solute center of geometries
    topology::Chargegroup_Iterator
      cg1 =   topo.chargegroup_begin();
    
    unsigned int i = 0, num_cg = topo.num_solute_chargegroups();
    
    for(i=0; i < num_cg; ++cg1, ++i){
      cg1.cog(pos, m_cg_cog(i));
    }

  } // chargegroup based cutoff

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// pairlist update
////////////////////////////////////////////////////////////////////////////////

void interaction::Standard_Pairlist_Algorithm::
update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       interaction::PairlistContainer & pairlist,
       unsigned int begin, unsigned int end,
       unsigned int stride)
{
  if (sim.param().pairlist.atomic_cutoff){
    // see standard_pairlist_algorithm_atomic.cc
    update_atomic(topo, conf, sim, pairlist, begin, end, stride);
  }
  else{
    update_cg(topo, conf, sim, pairlist, begin, end, stride);
  }
}

void interaction::Standard_Pairlist_Algorithm::update_cg(
  topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  interaction::PairlistContainer & pairlist,
  unsigned int begin, unsigned int end,
  unsigned int stride) {
  SPLIT_BOUNDARY(_update_cg, topo, conf, sim,
    	          pairlist, begin, end, stride);
}

//template<typename t_interaction_spec>
template<math::boundary_enum b>
void interaction::Standard_Pairlist_Algorithm::_update_cg
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::PairlistContainer & pairlist,
 unsigned int begin, unsigned int end,
 unsigned int stride
 )
{

  DEBUG(7, "standard pairlist update");
  if (begin == 0) 
    timer().start("pairlist");
  
  math::Periodicity<b> periodicity(conf.current().box);
  // empty the pairlist
  pairlist.clear();

  // loop over the chargegroups
  const int num_solute_cg = topo.num_solute_chargegroups();
  const int num_cg = topo.num_chargegroups();
  DEBUG(8, "begin=" << begin << " stride=" << stride
	<< " num_cg=" << num_cg << "num_solute_cg=" << num_solute_cg);

  int cg1 = 0, cg2 = 0;
  math::Vec r;
  
  const simulation::qmmm_enum qmmm = sim.param().qmmm.qmmm;
  
  // solute -
  for(cg1 = begin; cg1 < num_solute_cg; cg1+=stride){
    
    DEBUG(10, "cg1 = " << cg1);
    
    // If cg is QM
    if (Pairlist_Algorithm::qm_excluded(topo, qmmm, topo.chargegroup(cg1))) {
      DEBUG(9, "Skipping all for cg " << cg1);
      DEBUG(9, " - atoms " << topo.chargegroup(cg1) << "-" << topo.chargegroup(cg1+1)-1);
      continue;
    }

    if (!qmmm || !topo.is_qm( topo.chargegroup(cg1) )) { // skip QM chargegroups
      for(int a1 = topo.chargegroup(cg1),
      a_to = topo.chargegroup(cg1+1);
    a1 < a_to; ++a1){
        for(int a2 = a1+1; a2 < a_to; ++a2){
    
    // check it is not excluded
    if (excluded_solute_pair(topo, a1, a2))
      continue;
    
    assert(int(pairlist.size()) > a1);
    pairlist.solute_short[a1].push_back(a2);
        }
      }
    }
    else {
      DEBUG(9, "Skipping cg " << cg1 << " innerloop");
      DEBUG(9, " - atoms " << topo.chargegroup(cg1) << "-" << topo.chargegroup(cg1+1)-1);
    }
    
    // solute - solute
    DEBUG(10, "solute - solute");
    
    for(cg2 = cg1+1; cg2 < num_solute_cg; ++cg2){

      DEBUG(10, "cg2 = " << cg2);
      // If cg is QM
      if (Pairlist_Algorithm::qm_excluded(
            topo, qmmm, topo.chargegroup(cg1), topo.chargegroup(cg2))) 
        {
        DEBUG(9, "Skipping cgs " << cg1 << " and " << cg2);
        DEBUG(9, " - atoms " << topo.chargegroup(cg1) << "-" << topo.chargegroup(cg1+1)-1);
        DEBUG(9, " - atoms " << topo.chargegroup(cg2) << "-" << topo.chargegroup(cg2+1)-1);
        continue;
      }
      assert(m_cg_cog.size() > unsigned(cg1) &&
	     m_cg_cog.size() > unsigned(cg2));
      DEBUG(10, "ni cog1"<< math::v2s(m_cg_cog(cg1)));
      DEBUG(10, "ni cog2"<< math::v2s(m_cg_cog(cg2)));
      periodicity.nearest_image(m_cg_cog(cg1), m_cg_cog(cg2), r);
      DEBUG(10, "ni r"<< math::v2s(r));
      // the distance
      const double d = math::abs2(r);
      
      if (d > m_cutoff_long_2){        // OUTSIDE
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
      
	for(int a1 = topo.chargegroup(cg1),
	      a1_to = topo.chargegroup(cg1+1);
	    a1 != a1_to; ++a1){
	  for(int a2 = topo.chargegroup(cg2),
		a2_to = topo.chargegroup(cg2+1);
	      a2 != a2_to; ++a2){

	    pairlist.solute_long[a1].push_back(a2);

	  } // loop over atom of cg2
	} // loop over atom of cg1
	
	continue;
      } // longrange

      // SHORTRANGE : at least the second cg is solvent => no exclusions
      DEBUG(15, "a1=" << topo.chargegroup(cg1) << " a1_to=" << topo.chargegroup(cg1+1));
      DEBUG(15, "a2=" << topo.chargegroup(cg2) << " a2_to=" << topo.chargegroup(cg2+1));

      for(int a1 = topo.chargegroup(cg1),
	    a1_to = topo.chargegroup(cg1+1);
	  a1 != a1_to; ++a1){

	assert(a1 < int(topo.num_solute_atoms()));
	
	for(int a2 = topo.chargegroup(cg2),
	      a2_to = topo.chargegroup(cg2+1);
	    a2 != a2_to; ++a2){
	  
	  assert(a2 < int(topo.num_solute_atoms()));
	  DEBUG(16, "checking excl " << a1 << " - " << a2);
	  
	  if (excluded_solute_pair(topo, a1, a2)){
	    DEBUG(16, "->excluded (continue)");
	    continue;
	  }
	  
	  assert(int(pairlist.size()) > a1);
	  DEBUG(16, "push back a1=" << a1 << " in a pairlist " << pairlist.size());
	  pairlist.solute_short[a1].push_back(a2);
	  DEBUG(16, "a2=" << a2 << " done");
	  
	} // loop over atom of cg2
	DEBUG(16, "a2 for a1=" << a1 << " done");
      } // loop over atom of cg1
      
    }

    // solute - solvent
    DEBUG(10, "solute - solvent");
    
    for( ; cg2 < num_cg; ++cg2){

      DEBUG(10, "cg2 = " << cg2);// If cg is QM
      if (Pairlist_Algorithm::qm_excluded(
            topo, qmmm, topo.chargegroup(cg1), topo.chargegroup(cg2))) 
        {
        DEBUG(9, "Skipping cgs " << cg1 << " and " << cg2);
        DEBUG(9, " - atoms " << topo.chargegroup(cg1) << "-" << topo.chargegroup(cg1+1)-1);
        DEBUG(9, " - atoms " << topo.chargegroup(cg2) << "-" << topo.chargegroup(cg2+1)-1);
        continue;
      }
      
      assert(m_cg_cog.size() > unsigned(cg1));
    
      periodicity.nearest_image(m_cg_cog(cg1), 
				conf.current().pos(topo.chargegroup(cg2)),
				r);
    
      // the distance
      const double d = math::abs2(r);
      
      if (d > m_cutoff_long_2){        // OUTSIDE
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	
	for(int a1 = topo.chargegroup(cg1),
	      a1_to = topo.chargegroup(cg1+1);
	    a1 != a1_to; ++a1){
	  for(int a2 = topo.chargegroup(cg2),
		a2_to = topo.chargegroup(cg2+1);
	      a2 != a2_to; ++a2){
	    
	    // the interactions
	    pairlist.solute_long[a1].push_back(a2);
	    
	  } // loop over atom of cg2
	} // loop over atom of cg1
	
	continue;
      } // longrange
      
      // SHORTRANGE : at least the second cg is solvent => no exclusions
      for(int a1 = topo.chargegroup(cg1),
	    a1_to = topo.chargegroup(cg1+1);
	  a1 != a1_to; ++a1){
	for(int a2 = topo.chargegroup(cg2),
	      a2_to = topo.chargegroup(cg2+1);
	    a2 != a2_to; ++a2){
	  
	  pairlist.solute_short[a1].push_back(a2);
	  
	} // loop over atom of cg2
      } // loop over atom of cg1
    }
    
  } // cg1

  DEBUG(10, "solvent - solvent");

  // solvent - solvent
  const bool no_cuda = sim.param().innerloop.method != simulation::sla_cuda;
  if (no_cuda && num_cg > num_solute_cg)
    _solvent_solvent(topo, conf, sim, pairlist, cg1, stride, periodicity);

  if (begin == 0)
    timer().stop("pairlist");
  DEBUG(7, "pairlist done");

}

template<math::boundary_enum b>
void interaction::Standard_Pairlist_Algorithm::_solvent_solvent
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::PairlistContainer & pairlist,
 int cg1,  int stride,
 math::Periodicity<b> const & periodicity
 )
{
  bool master = false;
  if (cg1 == int(topo.num_solute_chargegroups())) { // master
    timer().start("pairlist solvent-solvent");
    master = true;
  }
  
  const int num_cg = topo.num_chargegroups();
  math::Vec r;
  
  for( ; cg1 < num_cg; cg1+=stride){

    for(int cg2 = cg1+1; cg2 < num_cg; ++cg2){

      periodicity.nearest_image(conf.current().pos(topo.chargegroup(cg1)),
				conf.current().pos(topo.chargegroup(cg2)),
				r);
    
      // the distance
      const double d = math::abs2(r);
    
      if (d > m_cutoff_long_2){        // OUTSIDE
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	
	for(int a1 = topo.chargegroup(cg1),
	      a1_to = topo.chargegroup(cg1+1);
	    a1 != a1_to; ++a1){
	  for(int a2 = topo.chargegroup(cg2),
		a2_to = topo.chargegroup(cg2+1);
	      a2 != a2_to; ++a2){
            DEBUG(15, "solvent-solvent longrange: " << a1 << "-" << a2);
	    pairlist.solvent_long[a1].push_back(a2);
	    
	  } // loop over atom of cg2
	} // loop over atom of cg1
	
	continue;
      } // longrange
      
      // SHORTRANGE : at least the second cg is solvent => no exclusions
      for(int a1 = topo.chargegroup(cg1),
	    a1_to = topo.chargegroup(cg1+1);
	  a1 != a1_to; ++a1){
	for(int a2 = topo.chargegroup(cg2),
	      a2_to = topo.chargegroup(cg2+1);
	    a2 != a2_to; ++a2){
	  DEBUG(15, "solvent-solvent shortrange: " << a1 << "-" << a2);
	  pairlist.solvent_short[a1].push_back(a2);
	  
	} // loop over atom of cg2
      } // loop over atom of cg1
      
    }
    
  } // cg1

  if (master)
    timer().stop("pairlist solvent-solvent");
  
}

////////////////////////////////////////////////////////////////////////////////
// perturbation
////////////////////////////////////////////////////////////////////////////////


void interaction::Standard_Pairlist_Algorithm::
update_perturbed(topology::Topology & topo,
		 configuration::Configuration & conf,
		 simulation::Simulation & sim,
		 interaction::PairlistContainer & pairlist,
		 interaction::PairlistContainer & perturbed_pairlist,
		 unsigned int begin, unsigned int end,
		 unsigned int stride)
{
  if (sim.param().pairlist.atomic_cutoff){
    update_perturbed_atomic(topo, conf, sim, 
			    pairlist, perturbed_pairlist,
			    begin, end, stride);
  }
  else{
    SPLIT_BOUNDARY(_update_pert_cg,
			 topo, conf, sim,
			 pairlist, perturbed_pairlist, 
			 begin, end, stride);
  }
}

template<math::boundary_enum b>
void interaction::Standard_Pairlist_Algorithm::
_update_pert_cg(topology::Topology & topo,
		configuration::Configuration & conf,
		simulation::Simulation & sim,
		interaction::PairlistContainer & pairlist,
		interaction::PairlistContainer & perturbed_pairlist,
		unsigned int begin, unsigned int end,
		unsigned int stride)
{
  DEBUG(7, "standard pairlist update");
  if (begin == 0) // master
    timer().start("perturbed pairlist");
  
  // create the innerloops
  math::Periodicity<b> periodicity(conf.current().box);

  // check whether we do scaling && scaling only
  bool scaled_only = (sim.param().perturbation.scaling && sim.param().perturbation.scaled_only);
  // empty the pairlist
  assert(pairlist.size() == topo.num_atoms());
  assert(perturbed_pairlist.size() == topo.num_atoms());

  pairlist.clear();
  perturbed_pairlist.clear();

  DEBUG(7, "pairlist cleared");

  // loop over the chargegroups
  const int num_solute_cg = topo.num_solute_chargegroups();
  const int num_cg = topo.num_chargegroups();
  DEBUG(8, "begin=" << begin << " stride=" << stride
	<< " num_cg=" << num_cg << "num_solute_cg=" << num_solute_cg);

  int cg1 = 0, cg2 = 0;
  math::Vec r;
  
  // solute -
  for(cg1 = begin; cg1 < num_solute_cg; cg1+=stride){
    
    for(int a1 = topo.chargegroup(cg1),
	  a_to = topo.chargegroup(cg1+1);
	a1 < a_to; ++a1){
      for(int a2 = a1+1; a2 < a_to; ++a2){
	
	// check it is not excluded
	if (excluded_solute_pair(topo, a1, a2))
	  continue;

	if (insert_pair(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
						a1, a2, scaled_only))
	  ;
	else if (insert_pair(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
						     a2, a1, scaled_only))
	  ;
	else
	  pairlist.solute_short[a1].push_back(a2);
	
      }
    }
    
    // solute - solute
    for(cg2 = cg1+1; cg2 < num_solute_cg; ++cg2){

      assert(m_cg_cog.size() > unsigned(cg1) &&
	     m_cg_cog.size() > unsigned(cg2));
    
      periodicity.nearest_image(m_cg_cog(cg1), m_cg_cog(cg2), r);

      // the distance
      const double d = math::abs2(r);
    
      if (d > m_cutoff_long_2){        // OUTSIDE
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE
      
	for(int a1 = topo.chargegroup(cg1),
	      a1_to = topo.chargegroup(cg1+1);
	    a1 != a1_to; ++a1){
	  for(int a2 = topo.chargegroup(cg2),
		a2_to = topo.chargegroup(cg2+1);
	      a2 != a2_to; ++a2){

	    if (insert_pair(topo, pairlist.solute_long, perturbed_pairlist.solute_long,
						a1, a2, scaled_only))
              {}
	    else if (insert_pair(topo, pairlist.solute_long, perturbed_pairlist.solute_long,
						a2, a1, scaled_only))
	      {}
	    else
	      pairlist.solute_long[a1].push_back(a2);

	  } // loop over atom of cg2
	} // loop over atom of cg1
	
	continue;
      } // longrange

      // SHORTRANGE : at least the second cg is solvent => no exclusions
	for(int a1 = topo.chargegroup(cg1),
	      a1_to = topo.chargegroup(cg1+1);
	    a1 != a1_to; ++a1){
	  for(int a2 = topo.chargegroup(cg2),
		a2_to = topo.chargegroup(cg2+1);
	      a2 != a2_to; ++a2){

	    if (excluded_solute_pair(topo, a1, a2))
	      continue;

	    if (insert_pair(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
						    a1, a2, scaled_only))
	      {}
	    else if (insert_pair(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
							 a2, a1, scaled_only))
	      {}
	    else
	      pairlist.solute_short[a1].push_back(a2);

	  } // loop over atom of cg2
	} // loop over atom of cg1
      
    }

    // solute - solvent
    for( ; cg2 < num_cg; ++cg2){

      assert(m_cg_cog.size() > unsigned(cg1));
    
      periodicity.nearest_image(m_cg_cog(cg1), 
				conf.current().pos(topo.chargegroup(cg2)),
				r);
    
      // the distance
      const double d = math::abs2(r);
    
      if (d > m_cutoff_long_2){        // OUTSIDE
	continue;
      }
  
      if (d > m_cutoff_short_2){       // LONGRANGE
	
	for(int a1 = topo.chargegroup(cg1),
	      a1_to = topo.chargegroup(cg1+1);
	    a1 != a1_to; ++a1){
	  for(int a2 = topo.chargegroup(cg2),
		a2_to = topo.chargegroup(cg2+1);
	      a2 != a2_to; ++a2){
	    
	    if (insert_pair(topo, pairlist.solute_long, perturbed_pairlist.solute_long,
						    a1, a2, scaled_only))
	      {}
	    else
	      pairlist.solute_long[a1].push_back(a2);
	    
	  } // loop over atom of cg2
	} // loop over atom of cg1
	
	continue;
      } // longrange
      
      // SHORTRANGE : at least the second cg is solvent => no exclusions
      for(int a1 = topo.chargegroup(cg1),
	    a1_to = topo.chargegroup(cg1+1);
	  a1 != a1_to; ++a1){
	for(int a2 = topo.chargegroup(cg2),
	      a2_to = topo.chargegroup(cg2+1);
	    a2 != a2_to; ++a2){
	  
	  if (insert_pair(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
						  a1, a2, scaled_only))
	    ;
	  else
	    pairlist.solute_short[a1].push_back(a2);
	  
	} // loop over atom of cg2
      } // loop over atom of cg1
      
    }
    
  } // cg1

  // solvent - solvent
  
  const bool no_cuda = sim.param().innerloop.method != simulation::sla_cuda; //--martina
  if (no_cuda && num_cg > num_solute_cg)
    _solvent_solvent(topo, conf, sim, pairlist, cg1, stride, periodicity); //only create solvent-solvent pairlist if we dont have cuda.
  
  /*_solvent_solvent(topo, conf, sim,
                   pairlist,  
                   cg1, stride, periodicity);
  */
  if (begin == 0) // master
    timer().stop("perturbed pairlist");
  DEBUG(7, "pairlist done");

}

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

bool interaction::Standard_Pairlist_Algorithm
::excluded_solute_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j)
{
  assert(i<j);
  return topo.all_exclusion(i).is_excluded(j);
}

