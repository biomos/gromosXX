/**
 * @file nonbonded_interaction.tcc
 * template methods of Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
template<typename t_interaction_spec>
inline
interaction::Nonbonded_Outerloop<t_interaction_spec>
::Nonbonded_Outerloop(Nonbonded_Parameter &nbp)
  : Nonbonded_Innerloop<t_interaction_spec>(nbp)
{
}

//==================================================
// interaction loops
//==================================================

/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop<t_interaction_spec>
::lj_crf_outerloop(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim, 
		   std::vector<std::vector<size_t> > const & pairlist,
		   Storage & storage)
{  
  DEBUG(7, "\tcalculate interactions");  

  Periodicity_type periodicity(conf.current().box);

  /*
    variables for a OMP parallelizable loop.
    outer index has to be integer...
  */
  std::vector<size_t>::const_iterator j_it, j_to;
  int i;
  int size_i = pairlist.size();

  //**************************
  // the Bekker implementation
  //**************************
  if (t_interaction_spec::do_bekker){

    periodicity.recalc_shift_vectors();

    int pc;
    size_t j;
    // translate the atom j
    DEBUG(9, "nonbonded_interaction: grid based pairlist");

    for(i=0; i < size_i; ++i){

      for(j_it = pairlist[i].begin(),
	    j_to = pairlist[i].end();
	  j_it != j_to;
	  ++j_it){
      
	pc = (*j_it >> 26);
	j = (*j_it & 67108863);
      
	DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << j
	      << " pc " << pc);
      
	lj_crf_innerloop(topo, conf, i, j, storage, periodicity, pc);
      }
      
    }

  }
  //*************************
  // standard implementation
  //*************************
  else{ // no grid based pairlist

    DEBUG(9, "nonbonded_interaction: no grid based pairlist");
    for(i=0; i < size_i; ++i){
    
      for(j_it = pairlist[i].begin(),
	    j_to = pairlist[i].end();
	  j_it != j_to;
	  ++j_it){

	DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << *j_it);
	// printf("nb pair %d - %d\n", i, *j_it);
	
	// shortrange, therefore store in simulation.system()
	lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);

	// storage.energies.bond_energy[0] += *j_it;

      }
      
    }
  }
}


/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop<t_interaction_spec>
::one_four_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage)
{
  DEBUG(7, "\tcalculate 1,4-interactions");

  Periodicity_type periodicity(conf.current().box);

  std::set<int>::const_iterator it, to;
  size_t const num_solute_atoms = topo.num_solute_atoms();
  size_t i;
  
  for(i=0; i < num_solute_atoms; ++i){
    it = topo.one_four_pair(i).begin();
    to = topo.one_four_pair(i).end();
    
    for( ; it != to; ++it){

      one_four_interaction_innerloop(topo, conf, i, *it, periodicity);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop<t_interaction_spec>
::RF_excluded_outerloop(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			Storage & storage)
{
  
  DEBUG(7, "\tcalculate RF excluded interactions");

  Periodicity_type periodicity(conf.current().box);
  
  int i, num_solute_atoms = topo.num_solute_atoms();
  
  for(i=0; i<num_solute_atoms; ++i){
    
    RF_excluded_interaction_innerloop(topo, conf, i, periodicity);
    
  } // loop over solute atoms

  // Solvent
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();
  cg_it += topo.num_solute_chargegroups();

  for(; cg_it != cg_to; ++cg_it){

    RF_solvent_interaction_innerloop(topo, conf, cg_it, periodicity);
    ++cg_it;

  } // loop over solvent charge groups

}  

