/**
 *  standard_pairlist_algorithm.tcc
 * create an atomic pairlist with a
 * chargegroup or an atom based cut-off criterion.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include <util/debug.h>

template<typename t_nonbonded_spec>
inline
interaction::Standard_Pairlist_Algorithm<t_nonbonded_spec>::
Standard_Pairlist_Algorithm()
  : interaction::Pairlist_Algorithm<t_nonbonded_spec>(),
    t_nonbonded_spec::exclusion_filter_type(),
    t_nonbonded_spec::range_filter_type()
{
}

template<typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_nonbonded_spec>::
update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim, 
       t_nonbonded_interaction & nonbonded_interaction)
{
  DEBUG(7, "pairlist update");

  Periodicity_type periodicity(conf.current().box);
   
  // empty the pairlist
  nonbonded_interaction.pairlist().clear();
  nonbonded_interaction.pairlist().resize(topo.num_atoms());

  if(t_nonbonded_spec::do_perturbation){
    // and the perturbed pairlist
    nonbonded_interaction.perturbed_pairlist().clear();
    nonbonded_interaction.perturbed_pairlist().resize(topo.num_atoms());
  }
  
  DEBUG(7, "pairlist(s) resized");
  
  // prepare the range filter (cutoff)
  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);
  
  if (!t_nonbonded_spec::do_atomic_cutoff){
    // prepare the range filter (center of geometries)    
    prepare_cog(topo, conf, sim);
    DEBUG(7, "range filter prepared (cog)");

  }
  
  // loop over the chargegroups
  topology::Chargegroup_Iterator cg1 = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();
  
  /**
   * @todo rewrite into two loops (one solute, one solvent)
   * or maybe not (for easier parallelization)
   */
  const int num_cg = topo.num_chargegroups();
  int cg1_index;
  
  for(cg1_index=0; cg1_index < num_cg; ++cg1_index) {
    // add intra cg (if not solvent...)
    
    cg1 = topo.chargegroup_it(cg1_index);

    if (unsigned(**cg1) < topo.solute().num_atoms()){
      do_cg_interaction_intra(topo, conf, sim, nonbonded_interaction, cg1,
			      periodicity);
    }
  
    // inter chargegroup
    topology::Chargegroup_Iterator cg2(*cg1+1);
    for(int cg2_index = cg1_index + 1; cg2 != cg_to; ++cg2, ++cg2_index) {

      if (!t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (range_chargegroup_pair(topo, conf, sim, nonbonded_interaction,
				   cg1_index, cg2_index, cg1, cg2,
				   periodicity))
	  continue;
      }
      
      DEBUG(11, "\tshortrange!");
      
      // SHORTRANGE
      if (unsigned(**cg2) < topo.solute().num_atoms()){
	// exclusions! (because cg2 is not solvent)
	do_cg_interaction_excl(topo, conf, sim, 
			       nonbonded_interaction, cg1, cg2,
			       periodicity);
      }
      else{
	// no exclusions... (at least cg2 is solvent)
	do_cg_interaction(topo, conf, sim, 
			  nonbonded_interaction, cg1, cg2, periodicity);
      }
	
    } // inter cg (cg2)
  } // cg1
  
  DEBUG(7, "pairlist done");

}

/**
 * inter cg, no exclusion
 */
template<typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_nonbonded_spec>
::do_cg_interaction(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    t_nonbonded_interaction &nonbonded_interaction,
		    topology::Chargegroup_Iterator const &cg1,
		    topology::Chargegroup_Iterator const &cg2,
		    Periodicity_type const & periodicity,
		    int const pc)
{

  topology::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();

  DEBUG(11, "do_cg_interaction " << *a1);
    
  for( ; a1 != a1_to; ++a1){
    for(topology::Atom_Iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      if (t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (t_nonbonded_spec::do_bekker){
	  if (range_atom_pair(topo, conf, sim, 
			      nonbonded_interaction, *a1, *a2,
			      pc, periodicity))
	    continue;
	}
	else {
	  if (range_atom_pair(topo, conf, sim, 
			      nonbonded_interaction, *a1, *a2,
			      periodicity))
	    continue;

	}
      }
      
      if (t_nonbonded_spec::do_bekker)
	nonbonded_interaction.add_shortrange_pair(topo, conf, sim, 
						  *a1, *a2, pc);
      else
	nonbonded_interaction.add_shortrange_pair(topo, conf, sim, *a1, *a2);
    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}


template<typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_nonbonded_spec>
::do_cg_interaction_excl(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 t_nonbonded_interaction &nonbonded_interaction,
			 topology::Chargegroup_Iterator const & cg1,
			 topology::Chargegroup_Iterator const & cg2,
			 Periodicity_type const & periodicity,
			 int const pc)
{
  topology::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();

  DEBUG(11, "do_cg_interaction_excl " << *a1);
  
  for( ; a1 != a1_to; ++a1){
    for(topology::Atom_Iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      if (t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (t_nonbonded_spec::do_bekker){
	  if (range_atom_pair(topo, conf, sim, 
			      nonbonded_interaction, *a1, *a2,
			      pc, periodicity))
	    continue;
	}
	else{
	  if (range_atom_pair(topo, conf, sim, 
			      nonbonded_interaction, *a1, *a2,
			      periodicity))
	    continue;
	}
      }
      
      // check it is not excluded
      if (excluded_solute_pair(topo, conf, sim, *a1, *a2)){
	continue;
      }

      if (t_nonbonded_spec::do_bekker)
	nonbonded_interaction.add_shortrange_pair(topo, conf, sim, 
						  *a1, *a2, pc);
      else
	nonbonded_interaction.add_shortrange_pair(topo, conf, sim, 
						  *a1, *a2);
    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_nonbonded_spec>
::do_cg_interaction_inv_excl(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     t_nonbonded_interaction &nonbonded_interaction,
			     topology::Chargegroup_Iterator const & cg1,
			     topology::Chargegroup_Iterator const & cg2,
			     Periodicity_type const & periodicity, int const pc)
{
  topology::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();

  DEBUG(11, "do_cg_interaction_excl " << *a1);
  
  for( ; a1 != a1_to; ++a1){
    for(topology::Atom_Iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      if (t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (t_nonbonded_spec::do_bekker){
	  if (range_atom_pair(topo, conf, sim, nonbonded_interaction, 
			      *a1, *a2, pc, periodicity))
	    continue;
	}
	else{
	  if (range_atom_pair(topo, conf, sim, nonbonded_interaction, 
			      *a1, *a2, periodicity))
	    continue;
	}
      }

      // check it is not excluded
      if (inverse_excluded_solute_pair(topo, conf, sim, *a1, *a2)){
	continue;
      }

      if (t_nonbonded_spec::do_bekker)
	nonbonded_interaction.add_shortrange_pair(topo, conf, sim, 
						  *a1, *a2, pc);
      else
	nonbonded_interaction.add_shortrange_pair(topo, conf, sim, 
						  *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_nonbonded_spec>
::do_cg_interaction_intra(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  t_nonbonded_interaction & nonbonded_interaction,
			  topology::Chargegroup_Iterator const & cg1,
			  Periodicity_type const & periodicity, int const pc)
{
  topology::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();

  DEBUG(11, "do_cg_interaction_intra " << *a1);
  
  for( ; a1 != a1_to; ++a1){
    for(topology::Atom_Iterator
	  a2(*a1+1);
	a2 != a1_to; ++a2){

      // check it is not excluded
      if (excluded_solute_pair(topo, conf, sim, *a1, *a2))
	continue;

      if (t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (t_nonbonded_spec::do_bekker){
	  if (range_atom_pair(topo, conf, sim, nonbonded_interaction, 
			      *a1, *a2, pc, periodicity))
	    continue;
	}
	else{
	  if (range_atom_pair(topo, conf, sim, nonbonded_interaction, 
			      *a1, *a2, periodicity))
	    continue;
	}
      }

      if (t_nonbonded_spec::do_bekker)
	nonbonded_interaction.add_shortrange_pair(topo, conf, sim, 
						  *a1, *a2, pc);
      else
	nonbonded_interaction.add_shortrange_pair(topo, conf, sim, 
						  *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}
