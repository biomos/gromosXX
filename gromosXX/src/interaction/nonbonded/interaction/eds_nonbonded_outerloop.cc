/**
 * @file perturbed_nonbonded_outerloop.cc
 * (template) methods of Perturbed_Nonbonded_Outerloop.
 */
#ifdef XXMPI
#include <mpi.h>
#endif

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction_types.h"
#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/interaction/eds_nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_innerloop.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_innerloop.h"

#include "../../../interaction/nonbonded/interaction/eds_nonbonded_outerloop.h"

#include "../../../interaction/nonbonded/interaction_spec.h"

#include "../../../util/debug.h"
#include "../../../interaction/nonbonded/innerloop_template.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Eds_Nonbonded_Outerloop
::Eds_Nonbonded_Outerloop(Nonbonded_Parameter & nbp)
  : m_param(nbp)
{
}

//==================================================
// interaction loops
//==================================================

//==================================================
// the eds-perturbed interaction (outer) loops
//==================================================

/**
 * helper function to calculate perturbed forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
void interaction::Eds_Nonbonded_Outerloop
::eds_lj_crf_outerloop(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     Pairlist const & pairlist,
			     Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_eds_lj_crf_outerloop,
		       topo, conf, sim,
		       pairlist, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Eds_Nonbonded_Outerloop
::_eds_lj_crf_outerloop(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim,
			      Pairlist const & pairlist,
			      Storage & storage)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Eds_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);

  std::vector<unsigned int>::const_iterator j_it, j_to;
  unsigned int i = 0;
  unsigned int size_i = unsigned(topo.num_solute_atoms());
  // unsigned int size_i = unsigned(pairlist.size());

  DEBUG(6, "pert sr: " << size_i);

  for(i=0; i < size_i; ++i){
    
    for(j_it = pairlist[i].begin(),
	  j_to = pairlist[i].end();
	j_it != j_to;
	++j_it){
      
      DEBUG(10, "\tperturbed nonbonded_interaction: i "
	    << i << " j " << *j_it);
      
      innerloop.eds_lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
    
  }

  DEBUG(7, "end of function eds-perturbed nonbonded interaction");  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
void interaction::Eds_Nonbonded_Outerloop
::eds_one_four_outerloop(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_eds_one_four_outerloop,
		       topo, conf, sim, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Eds_Nonbonded_Outerloop
::_eds_one_four_outerloop(topology::Topology & topo,
				configuration::Configuration & conf,
				simulation::Simulation & sim,
				Storage & storage)
{
  DEBUG(7, "\tcalculate eds-perturbed 1,4-interactions");
  
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Eds_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);
  
  topology::excl_cont_t::value_type::const_iterator it, to;
  
  std::map<unsigned int, topology::EDS_Perturbed_Atom>::const_iterator 
    mit=topo.eds_perturbed_solute().atoms().begin(), 
    mto=topo.eds_perturbed_solute().atoms().end();
  
  for(; mit!=mto; ++mit){
    //int seq = mit->second.sequence_number();
    //it = topo.one_four_pair(seq).begin();
    //to = topo.one_four_pair(seq).end();
    
    it = mit->second.one_four_pair().begin();
    to = mit->second.one_four_pair().end();
    
    for( ; it != to; ++it){

      innerloop.eds_one_four_interaction_innerloop
	(topo, conf, mit->second.sequence_number(), *it, periodicity);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  


/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
void interaction::Eds_Nonbonded_Outerloop
::eds_RF_excluded_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_eds_RF_excluded_outerloop,
		       topo, conf, sim, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Eds_Nonbonded_Outerloop
::_eds_RF_excluded_outerloop(topology::Topology & topo,
				   configuration::Configuration & conf,
				   simulation::Simulation & sim,
				   Storage & storage)
{

  DEBUG(7, "\tcalculate eds-perturbed excluded RF interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Eds_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);

  std::map<unsigned int, topology::EDS_Perturbed_Atom>::const_iterator
    mit=topo.eds_perturbed_solute().atoms().begin(),
    mto=topo.eds_perturbed_solute().atoms().end();

  DEBUG(9, "\tSize of perturbed atoms " 
	<< unsigned(topo.eds_perturbed_solute().atoms().size()));
  
  for(; mit!=mto; ++mit){
    innerloop.eds_RF_excluded_interaction_innerloop(topo, conf, mit, periodicity);
  }
}



