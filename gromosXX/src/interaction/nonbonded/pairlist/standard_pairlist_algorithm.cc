/**
 * @file standard_pairlist_algorithm.cc
 * standard pairlist algorithm
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
#include <interaction/nonbonded/innerloop_template.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

interaction::Standard_Pairlist_Algorithm::
Standard_Pairlist_Algorithm()
  : interaction::Pairlist_Algorithm(),
		 m_solvent_solvent_timing(0.0),
		 m_spc_timing(0.0)
{
}

/**
 * put the chargegroups into the box
 */
template<math::boundary_enum b>
static void _prepare_cog(configuration::Configuration & conf,
			 topology::Topology & topo)
{
  DEBUG(10, "putting chargegroups into box");
  math::Periodicity<b> periodicity(conf.current().box);
  periodicity.put_chargegroups_into_box(conf, topo);
}

/**
 * calculate center of geometries
 */
void interaction::Standard_Pairlist_Algorithm::
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
      cg1 =   topo.chargegroup_begin(),
      cg_to = topo.chargegroup_end();
    
    unsigned int i, num_cg = topo.num_solute_chargegroups();
    
    for(i=0; i < num_cg; ++cg1, ++i){
      cg1.cog(pos, m_cg_cog(i));
    }

  } // chargegroup based cutoff
  
}

////////////////////////////////////////////////////////////////////////////////
// pairlist update
////////////////////////////////////////////////////////////////////////////////

void interaction::Standard_Pairlist_Algorithm::
update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       interaction::Storage & storage,
       interaction::Pairlist & pairlist,
       unsigned int begin, unsigned int end,
       unsigned int stride)
{
  if (sim.param().pairlist.atomic_cutoff){
    // see standard_pairlist_algorithm_atomic.cc
    update_atomic(topo, conf, sim, storage, pairlist, begin, end, stride);
  }
  else{
    update_cg(topo, conf, sim, storage, pairlist, begin, end, stride);
  }
}


void interaction::Standard_Pairlist_Algorithm::update_cg
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim, 
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 unsigned int begin, unsigned int end,
 unsigned int stride
 )
{
  SPLIT_INNERLOOP(_update_cg, topo, conf, sim, storage,
		  pairlist, begin, end, stride);
}

template<typename t_interaction_spec>
void interaction::Standard_Pairlist_Algorithm::_update_cg
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 unsigned int begin, unsigned int end,
 unsigned int stride
 )
{

  DEBUG(7, "standard pairlist update");
  const double update_start = util::now();
  
  // create the innerloop
  Nonbonded_Innerloop innerloop(*m_param);
  innerloop.init(sim);

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);

  // empty the pairlist
  for(unsigned int i=0; i<topo.num_atoms(); ++i)
    pairlist[i].clear();

  // loop over the chargegroups
  const int num_solute_cg = topo.num_solute_chargegroups();
  const int num_cg = topo.num_chargegroups();
  DEBUG(8, "begin=" << begin << " stride=" << stride
	<< " num_cg=" << num_cg << "num_solute_cg=" << num_solute_cg);

  int cg1, cg2;
  math::Vec r;
  
  // solute -
  for(cg1 = begin; cg1 < num_solute_cg; cg1+=stride){
    
    DEBUG(10, "cg1 = " << cg1);
    
    for(int a1 = topo.chargegroup(cg1),
	  a_to = topo.chargegroup(cg1+1);
	a1 < a_to; ++a1){
      for(int a2 = a1+1; a2 < a_to; ++a2){
	
	// check it is not excluded
	if (excluded_solute_pair(topo, a1, a2))
	  continue;
	
	assert(pairlist.size() > a1);
	pairlist[a1].push_back(a2);
	
      }
    }
    
    // solute - solute
    DEBUG(10, "solute - solute");
    
    for(cg2 = cg1+1; cg2 < num_solute_cg; ++cg2){

      DEBUG(10, "cg2 = " << cg2);
      
      assert(m_cg_cog.size() > unsigned(cg1) &&
	     m_cg_cog.size() > unsigned(cg2));
    
      periodicity.nearest_image(m_cg_cog(cg1), m_cg_cog(cg2), r);
    
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
	    innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);

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

	assert(a1 < topo.num_solute_atoms());
	
	for(int a2 = topo.chargegroup(cg2),
	      a2_to = topo.chargegroup(cg2+1);
	    a2 != a2_to; ++a2){
	  
	  assert(a2 < topo.num_solute_atoms());
	  DEBUG(16, "checking excl " << a1 << " - " << a2);
	  
	  if (excluded_solute_pair(topo, a1, a2)){
	    DEBUG(16, "->excluded (continue)");
	    continue;
	  }
	  
	  assert(pairlist.size() > a1);
	  DEBUG(16, "push back a1=" << a1 << " in a pairlist " << pairlist.size());
	  pairlist[a1].push_back(a2);
	  DEBUG(16, "a2=" << a2 << " done");
	  
	} // loop over atom of cg2
	DEBUG(16, "a2 for a1=" << a1 << " done");
      } // loop over atom of cg1
      
    }

    // solute - solvent
    DEBUG(10, "solute - solvent");
    
    for( ; cg2 < num_cg; ++cg2){

      DEBUG(10, "cg2 = " << cg2);
      
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
	    innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);
	    
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
	  
	  pairlist[a1].push_back(a2);
	  
	} // loop over atom of cg2
      } // loop over atom of cg1
    }
    
  } // cg1

  DEBUG(10, "solvent - solvent");

  // solvent - solvent
  if (sim.param().force.spc_loop)
    _spc_loop<t_interaction_spec>(topo, conf, sim, storage, pairlist, innerloop, cg1, stride, periodicity);    
  else
    _solvent_solvent<t_interaction_spec>(topo, conf, sim, storage, pairlist, innerloop, cg1, stride, periodicity);

  this->m_timing += util::now() - update_start;
  DEBUG(7, "pairlist done");

}

template<typename t_interaction_spec>
void interaction::Standard_Pairlist_Algorithm::_solvent_solvent
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 Nonbonded_Innerloop & innerloop,
 int cg1,  int stride,
 math::Periodicity<t_interaction_spec::boundary_type> const & periodicity
 )
{
  const double start = util::now();
  
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
	    
	    // the interactions
	    innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);
	    
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
	  
	  pairlist[a1].push_back(a2);
	  
	} // loop over atom of cg2
      } // loop over atom of cg1
      
    }
    
  } // cg1

  m_solvent_solvent_timing += util::now() - start;
  
}

template<typename t_interaction_spec>
void interaction::Standard_Pairlist_Algorithm::_spc_loop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 Nonbonded_Innerloop & innerloop,
 int cg1,  int stride,
 math::Periodicity<t_interaction_spec::boundary_type> const & periodicity
 )
{
  const double start = util::now();
  
  const int num_cg = topo.num_chargegroups();
  const int num_solute_cg = topo.num_solute_chargegroups();
  math::Vec r;

  int pl_index = topo.num_solute_atoms() + cg1 - topo.num_solute_chargegroups();
  
  DEBUG(7, "spc loop: solute cg=" << num_solute_cg << " cg=" << num_cg
	<< " cg1=" << cg1 << " pl index=" << pl_index << " stride=" << stride);

  for( ; cg1 < num_cg; cg1+=stride, pl_index += stride){

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
	
	innerloop.spc_innerloop<t_interaction_spec>(topo, conf, 
						    cg1 - num_solute_cg,
						    cg2 - num_solute_cg,
						    storage, periodicity);
	
	continue;
      } // longrange
      
      // SHORTRANGE : at least the second cg is solvent => no exclusions
      pairlist[pl_index].push_back(cg2 - num_solute_cg);
      
    }
    
  } // cg1

  m_spc_timing += util::now() - start;

}


////////////////////////////////////////////////////////////////////////////////
// perturbation
////////////////////////////////////////////////////////////////////////////////


void interaction::Standard_Pairlist_Algorithm::
update_perturbed(topology::Topology & topo,
		 configuration::Configuration & conf,
		 simulation::Simulation & sim,
		 interaction::Storage & storage,
		 interaction::Pairlist & pairlist,
		 interaction::Pairlist & perturbed_pairlist,
		 unsigned int begin, unsigned int end,
		 unsigned int stride)
{
  if (sim.param().pairlist.atomic_cutoff){
    update_perturbed_atomic(topo, conf, sim, storage,
			    pairlist, perturbed_pairlist,
			    begin, end, stride);
  }
  else{
    SPLIT_PERT_INNERLOOP(_update_pert_cg,
			 topo, conf, sim, storage,
			 pairlist, perturbed_pairlist, 
			 begin, end, stride);
  }
}

template<typename t_interaction_spec, typename t_perturbation_details>
void interaction::Standard_Pairlist_Algorithm::
_update_pert_cg(topology::Topology & topo,
		configuration::Configuration & conf,
		simulation::Simulation & sim,
		interaction::Storage & storage,
		interaction::Pairlist & pairlist,
		interaction::Pairlist & perturbed_pairlist,
		unsigned int begin, unsigned int end,
		unsigned int stride)
{
  DEBUG(7, "standard pairlist update");
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

  // loop over the chargegroups
  const int num_solute_cg = topo.num_solute_chargegroups();
  const int num_cg = topo.num_chargegroups();
  DEBUG(8, "begin=" << begin << " stride=" << stride
	<< " num_cg=" << num_cg << "num_solute_cg=" << num_solute_cg);

  int cg1, cg2;
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

	if (insert_pair<t_perturbation_details>(topo, pairlist, perturbed_pairlist,
						a1, a2, sim.param().perturbation.scaled_only))
	  ;
	else if (insert_pair<t_perturbation_details>(topo, pairlist, perturbed_pairlist,
						     a2, a1, sim.param().perturbation.scaled_only))
	  ;
	else
	  pairlist[a1].push_back(a2);
	
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
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
      
	for(int a1 = topo.chargegroup(cg1),
	      a1_to = topo.chargegroup(cg1+1);
	    a1 != a1_to; ++a1){
	  for(int a2 = topo.chargegroup(cg2),
		a2_to = topo.chargegroup(cg2+1);
	      a2 != a2_to; ++a2){

	    if (calculate_pair<t_interaction_spec, t_perturbation_details>
		(topo, conf, storage, innerloop, perturbed_innerloop, a1, a2,
		 periodicity, sim.param().perturbation.scaled_only))
	      {}
	    else if (calculate_pair<t_interaction_spec, t_perturbation_details>
		(topo, conf, storage, innerloop, perturbed_innerloop, a2, a1,
		 periodicity, sim.param().perturbation.scaled_only))
	      {}
	    else
	      innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);

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

	    if (insert_pair<t_perturbation_details>(topo, pairlist, perturbed_pairlist,
						    a1, a2, sim.param().perturbation.scaled_only))
	      {}
	    else if (insert_pair<t_perturbation_details>(topo, pairlist, perturbed_pairlist,
							 a2, a1, sim.param().perturbation.scaled_only))
	      {}
	    else
	      pairlist[a1].push_back(a2);

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
  
      if (d > m_cutoff_short_2){       // LONGRANGE: calculate interactions!
	
	for(int a1 = topo.chargegroup(cg1),
	      a1_to = topo.chargegroup(cg1+1);
	    a1 != a1_to; ++a1){
	  for(int a2 = topo.chargegroup(cg2),
		a2_to = topo.chargegroup(cg2+1);
	      a2 != a2_to; ++a2){
	    
	    // the interactions
	    if (calculate_pair<t_interaction_spec, t_perturbation_details>
		(topo, conf, storage, innerloop, perturbed_innerloop, a1, a2,
		 periodicity, sim.param().perturbation.scaled_only))
	      {}
	    else
	      innerloop.lj_crf_innerloop<t_interaction_spec>(topo, conf, a1, a2, storage, periodicity);
	    
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
	  
	  if (insert_pair<t_perturbation_details>(topo, pairlist, perturbed_pairlist,
						  a1, a2, sim.param().perturbation.scaled_only))
	    ;
	  else
	    pairlist[a1].push_back(a2);
	  
	} // loop over atom of cg2
      } // loop over atom of cg1
      
    }
    
  } // cg1

  // solvent - solvent
  if (sim.param().force.spc_loop)
    _spc_loop<t_interaction_spec>(topo, conf, sim,
				  storage, pairlist, innerloop, 
				  cg1, stride, periodicity);    
  else
    _solvent_solvent<t_interaction_spec>(topo, conf, sim,
					 storage, pairlist, innerloop, 
					 cg1, stride, periodicity);

  this->m_timing += util::now() - update_start;
  DEBUG(7, "pairlist done");

}


template<typename t_perturbation_details>
inline bool interaction::Standard_Pairlist_Algorithm::insert_pair
(
 topology::Topology & topo,
 interaction::Pairlist & pairlist,
 interaction::Pairlist & perturbed_pairlist,
 int a1, int a2,
 bool scaled_only
 )
{
  if (topo.is_perturbed(a1)){
    if (t_perturbation_details::do_scaling && scaled_only){
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


template<typename t_interaction_spec, typename t_perturbation_details>
inline bool interaction::Standard_Pairlist_Algorithm::calculate_pair
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 interaction::Storage & storage,
 Nonbonded_Innerloop & innerloop,
 Perturbed_Nonbonded_Innerloop
 <t_interaction_spec, t_perturbation_details> & perturbed_innerloop,
 int a1, int a2,
 math::Periodicity<t_interaction_spec::boundary_type> const & periodicity,
 bool scaled_only
 )
{

  if (topo.is_perturbed(a1)){
    if (t_perturbation_details::do_scaling && scaled_only){
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

    return true;
  }
  return false;
}

bool interaction::Standard_Pairlist_Algorithm
::excluded_solute_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j)
{
  assert(i<j);
  
  std::set<int>::const_reverse_iterator
    e = topo.all_exclusion(i).rbegin(),
    e_to = topo.all_exclusion(i).rend();

  for( ; e != e_to; ++e){
    if (j > unsigned(*e)) break;
    if (j == unsigned(*e)){
      DEBUG(11, "\texcluded");
      return true;
    }
      
  }
  DEBUG(12, "\tnot excluded");
  return false;
}
