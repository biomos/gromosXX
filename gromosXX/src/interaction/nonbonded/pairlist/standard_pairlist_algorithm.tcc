/**
 * @file standard_pairlist_algorithm.tcc
 * create an atomic pairlist with a
 * chargegroup or an atom based cut-off criterion.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

template<typename t_interaction_spec, typename t_perturbation_spec>
inline
interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>::
Standard_Pairlist_Algorithm()
  : interaction::Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>()
{
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>::
prepare(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{

  // prepare the range filter (cutoff and chargegroup center of geometries)
  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);
  
  if (!t_interaction_spec::do_atomic_cutoff){

    this->prepare_cog(topo, conf, sim);
    DEBUG(7, "range filter prepared (cog)");

  }
  
}


template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>::
update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim, 
       Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
       unsigned int begin, unsigned int end, unsigned int stride)
{
  DEBUG(7, "standard pairlist update");

  const double update_start = util::now();

  Periodicity_type periodicity(conf.current().box);
   
  // empty the pairlist
  // to be general, we anyway empty the COMPLETE pairlist
  // (changing begin, end, stride during simulation possible)
  for(size_t i=0; i<topo.num_atoms(); ++i)
    nbs.pairlist()[i].clear();

  if(t_perturbation_spec::do_perturbation){
    // and the perturbed pairlist
    for(size_t i=0; i<topo.num_atoms(); ++i)
      nbs.perturbed_pairlist()[i].clear();

  }
  
  DEBUG(7, "pairlist(s) resized");
  
  
  // loop over the chargegroups
  
  /**
   * @todo rewrite into two loops (one solute, one solvent)
   * or maybe not (for easier parallelization)
   */
  const int num_cg = topo.num_chargegroups();
  const int num_solute_cg = topo.num_solute_chargegroups();
  int cg1_index, cg1_to;

  topology::Chargegroup_Iterator cg1;

  cg1_index = begin;
  cg1_to = num_cg;
  for( ; cg1_index < num_solute_cg; cg1_index+=stride) {

    cg1 = topo.chargegroup_it(cg1_index);

    // intra chargegroup
    do_cg_interaction_intra(topo, conf, sim, nbs, 
			    cg1, periodicity);

    // inter chargegroup
    do_cg1_loop(topo, conf, sim, nbs, cg1, 
		cg1_index, num_solute_cg, num_cg,
		periodicity);
    
  } // cg1

  for( ; cg1_index < num_cg; cg1_index+=stride) {

    // solvent
    cg1 = topo.chargegroup_it(cg1_index);

    
    do_cg1_loop(topo, conf, sim, nbs, cg1, 
		cg1_index, num_solute_cg, num_cg,
		periodicity);
    
  } // cg1

  this->m_timing += util::now() - update_start;
  
  DEBUG(7, "pairlist done");

}

/**
 * loop over chargegroup 1
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>
::do_cg1_loop(topology::Topology & topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
	      topology::Chargegroup_Iterator const & cg1,
	      int cg1_index,
	      int const num_solute_cg,
	      int const num_cg,
	      Periodicity_type const & periodicity)
{
  
  // inter chargegroup
  topology::Chargegroup_Iterator cg2 = *cg1+1;

  // solute...
  int cg2_index;
  for(cg2_index = cg1_index + 1; cg2_index < num_solute_cg; ++cg2, ++cg2_index){
    
    if (!t_interaction_spec::do_atomic_cutoff){
      // filter out interactions based on chargegroup distances
      
      if (this->range_chargegroup_pair(topo, conf, sim, nbs,
				       cg1_index, cg2_index, cg1, cg2,
				       periodicity))
	
	continue;
      
    }
    
    // SHORTRANGE
    // exclusions! (because cg2 is not solvent)
    do_cg_interaction_excl(topo, conf, sim, nbs,
			   cg1, cg2, periodicity);
    
  } // inter cg (cg2 solute)
  // solvent...
  for(; cg2_index < num_cg; ++cg2, ++cg2_index) {
    
    if (!t_interaction_spec::do_atomic_cutoff){
      // filter out interactions based on chargegroup distances
      
      if (this->range_chargegroup_pair(topo, conf, sim, nbs,
				       cg1_index, cg2_index, cg1, cg2,
				       periodicity))
	continue;
      
    }
    
    // SHORTRANGE
    do_cg_interaction(topo, conf, sim, nbs,
		      cg1, cg2, periodicity);
    
  } // inter cg (cg2 solvent)
  
}

/**
 * inter cg, no exclusion
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>
::do_cg_interaction(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
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

      if (t_interaction_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (t_interaction_spec::do_bekker){
	  if (this->range_atom_pair(topo, conf, sim, 
				    nbs, *a1, *a2,
				    pc, periodicity))
	    continue;
	}
	else {
	  if (this->range_atom_pair(topo, conf, sim, 
				    nbs, *a1, *a2,
				    periodicity))
	    continue;

	}
      }
      
      if (t_interaction_spec::do_bekker)
	nbs.add_shortrange_pair(topo, conf, sim, *a1, *a2, pc);
      else
	nbs.add_shortrange_pair(topo, conf, sim, *a1, *a2);
    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}


template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>
::do_cg_interaction_excl(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
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

      if (t_interaction_spec::do_atomic_cutoff){
	// filter out interactions based on atom distances
	if (t_interaction_spec::do_bekker){
	  if (this->range_atom_pair(topo, conf, sim, nbs, *a1, *a2,
				    pc, periodicity))
	    continue;
	}
	else{
	  if (this->range_atom_pair(topo, conf, sim, nbs, *a1, *a2,
				    periodicity))
	    continue;
	}
      }
      
      // check it is not excluded
      if (this->excluded_solute_pair(topo, conf, sim, *a1, *a2)){
	continue;
      }

      if (t_interaction_spec::do_bekker)
	nbs.add_shortrange_pair(topo, conf, sim, *a1, *a2, pc);
      else
	nbs.add_shortrange_pair(topo, conf, sim, *a1, *a2);
    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>
::do_cg_interaction_inv_excl(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
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

      if (t_interaction_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (t_interaction_spec::do_bekker){
	  if (this->range_atom_pair(topo, conf, sim, nbs, 
				    *a1, *a2, pc, periodicity))
	    continue;
	}
	else{
	  if (this->range_atom_pair(topo, conf, sim, nbs, 
				    *a1, *a2, periodicity))
	    continue;
	}
      }

      // check it is not excluded
      if (this->inverse_excluded_solute_pair(topo, conf, sim, *a1, *a2)){
	continue;
      }

      if (t_interaction_spec::do_bekker)
	nbs.add_shortrange_pair(topo, conf, sim, *a1, *a2, pc);
      else
	nbs.add_shortrange_pair(topo, conf, sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>
::do_cg_interaction_intra(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
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
      if (this->excluded_solute_pair(topo, conf, sim, *a1, *a2))
	continue;

      if (t_interaction_spec::do_atomic_cutoff){
	// filter out interactions based on chargegroup distances
	if (t_interaction_spec::do_bekker){
	  if (this->range_atom_pair(topo, conf, sim, nbs,
				    *a1, *a2, pc, periodicity))
	    continue;
	}
	else{
	  if (this->range_atom_pair(topo, conf, sim, nbs,
				    *a1, *a2, periodicity))
	    continue;
	}
      }

      if (t_interaction_spec::do_bekker)
	nbs.add_shortrange_pair(topo, conf, sim, *a1, *a2, pc);
      else
	nbs.add_shortrange_pair(topo, conf, sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}
