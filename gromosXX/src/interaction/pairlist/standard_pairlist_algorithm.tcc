/**
 *  standard_pairlist_algorithm.tcc
 * create an atomic pairlist with a
 * chargegroup or an atom based cut-off criterion.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include "../../debug.h"

template<typename t_simulation, typename t_nonbonded_spec>
inline
interaction::Standard_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>::
Standard_Pairlist_Algorithm()
  : interaction::Pairlist_Algorithm<t_simulation, t_nonbonded_spec>(),
    t_nonbonded_spec::exclusion_filter_type(),
    t_nonbonded_spec::range_filter_type()
{
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>::
update(t_simulation &sim, t_nonbonded_interaction & nonbonded_interaction)
{
  DEBUG(7, "pairlist update");
   
  // empty the pairlist
  nonbonded_interaction.pairlist().clear();
  nonbonded_interaction.pairlist().resize(sim.topology().num_atoms());

  if(t_nonbonded_spec::do_perturbation){
    // and the perturbed pairlist
    nonbonded_interaction.perturbed_pairlist().clear();
    nonbonded_interaction.perturbed_pairlist().resize(sim.topology().num_atoms());
  }
  
  DEBUG(7, "pairlist(s) resized");
  
  // prepare the range filter (cutoff)
  set_cutoff(sim.nonbonded().cutoff_short(), sim.nonbonded().cutoff_long());
  
  if (!t_nonbonded_spec::do_atomic_cutoff){
    // prepare the range filter (center of geometries)    
    prepare_cog(sim);
    DEBUG(7, "range filter prepared (cog)");

  }
  
  // loop over the chargegroups
  simulation::chargegroup_iterator cg1 = sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();
  
  /**
   * @TODO rewrite into two loops (one solute, one solvent)
   */
  for(int cg1_index=0; cg1 != cg_to; ++cg1, ++cg1_index) {
    // add intra cg (if not solvent...)
    if (unsigned(**cg1) < sim.topology().solute().num_atoms()){
      do_cg_interaction_intra(sim, nonbonded_interaction, cg1);
    }
  
    // inter chargegroup
    simulation::chargegroup_iterator cg2(*cg1+1);
    for(int cg2_index = cg1_index + 1; cg2 != cg_to; ++cg2, ++cg2_index) {

      if (!t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (range_chargegroup_pair(sim, nonbonded_interaction,
				   cg1_index, cg2_index, cg1, cg2))
	  continue;
      }
      
      DEBUG(11, "\tshortrange!");
      
      // SHORTRANGE
      if (unsigned(**cg2) < sim.topology().solute().num_atoms()){
	// exclusions! (because cg2 is not solvent)
	do_cg_interaction_excl(sim, nonbonded_interaction, cg1, cg2);
      }
      else{
	// no exclusions... (at least cg2 is solvent)
	do_cg_interaction(sim, nonbonded_interaction, cg1, cg2);
      }
	
    } // inter cg (cg2)
  } // cg1
  
  DEBUG(7, "pairlist done");

}

/**
 * inter cg, no exclusion
 */
template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>
::do_cg_interaction(t_simulation & sim,
		    t_nonbonded_interaction &nonbonded_interaction,
		    simulation::chargegroup_iterator const &cg1,
		    simulation::chargegroup_iterator const &cg2)
{

  simulation::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();

  DEBUG(11, "do_cg_interaction " << *a1);
    
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      if (t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (range_atom_pair(sim, nonbonded_interaction, *a1, *a2))
	  continue;
      }

      nonbonded_interaction.add_shortrange_pair(sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}


template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>
::do_cg_interaction_excl(t_simulation & sim,
			 t_nonbonded_interaction &nonbonded_interaction,
			 simulation::chargegroup_iterator const & cg1,
			 simulation::chargegroup_iterator const & cg2)
{
  simulation::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();

  DEBUG(11, "do_cg_interaction_excl " << *a1);
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      if (t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (range_atom_pair(sim, nonbonded_interaction, *a1, *a2))
	  continue;
      }

      // check it is not excluded
      if (excluded_solute_pair(sim, *a1, *a2)){
	continue;
      }

      nonbonded_interaction.add_shortrange_pair(sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>
::do_cg_interaction_inv_excl(t_simulation & sim,
			     t_nonbonded_interaction &nonbonded_interaction,
			     simulation::chargegroup_iterator const & cg1,
			     simulation::chargegroup_iterator const & cg2)
{
  simulation::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();

  DEBUG(11, "do_cg_interaction_excl " << *a1);
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      if (t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (range_atom_pair(sim, nonbonded_interaction, *a1, *a2))
	  continue;
      }

      // check it is not excluded
      if (inverse_excluded_solute_pair(sim, *a1, *a2)){
	continue;
      }

      nonbonded_interaction.add_shortrange_pair(sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Standard_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>
::do_cg_interaction_intra(t_simulation & sim,
			  t_nonbonded_interaction & nonbonded_interaction,
			  simulation::chargegroup_iterator const & cg1)
{
  simulation::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();

  DEBUG(11, "do_cg_interaction_intra " << *a1);
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2(*a1+1);
	a2 != a1_to; ++a2){

      // check it is not excluded
      if (excluded_solute_pair(sim, *a1, *a2))
	continue;

      if (t_nonbonded_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (range_atom_pair(sim, nonbonded_interaction, *a1, *a2))
	  continue;
      }

      nonbonded_interaction.add_shortrange_pair(sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}
