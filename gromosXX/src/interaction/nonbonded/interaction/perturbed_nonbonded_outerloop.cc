/**
 * @file perturbed_nonbonded_outerloop.cc
 * template methods of Perturbed_Nonbonded_Outerloop.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
template<typename t_interaction_spec, typename t_perturbation_details>
inline
interaction::Perturbed_Nonbonded_Outerloop<t_interaction_spec,  
					   t_perturbation_details>
::Perturbed_Nonbonded_Outerloop(Nonbonded_Parameter & nbp)
  : Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details>(nbp)
{
}

//==================================================
// interaction loops
//==================================================

//==================================================
// the perturbed interaction (outer) loops
//==================================================

/**
 * helper function to calculate perturbed forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_interaction_spec, typename  t_perturbation_details>
inline void interaction::Perturbed_Nonbonded_Outerloop<
  t_interaction_spec,  t_perturbation_details>
::perturbed_lj_crf_outerloop(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     std::vector<std::vector<unsigned int> > & pl,
			     Storage & storage)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  Periodicity_type periodicity(conf.current().box);

  std::vector<unsigned int>::const_iterator j_it, j_to;
  unsigned int i;
  unsigned int size_i = unsigned(pl.size());

  if (t_interaction_spec::do_bekker){

    periodicity.recalc_shift_vectors();
    
    int pc;
    unsigned int j;

    // translate the atom j

    for(i=0; i < size_i; ++i){
      
      for(j_it = pl[i].begin(),
	    j_to = pl[i].end();
	  j_it != j_to;
	  ++j_it){
      
	pc = (*j_it >> 26);
	j = (*j_it & 67108863);
      
	DEBUG(10, "\tperturbed nonbonded_interaction: i " << i << " j " << j
	      << " pc " << pc);
      
	perturbed_lj_crf_innerloop(topo, conf, i, j, storage, periodicity, pc);
      }
    }
  }
  else{

    DEBUG(9, "perturbed nonbonded_interaction: no grid based pairlist");
    for(i=0; i < size_i; ++i){
    
      for(j_it = pl[i].begin(),
	    j_to = pl[i].end();
	  j_it != j_to;
	  ++j_it){

	DEBUG(10, "\tperturbed nonbonded_interaction: i "
	      << i << " j " << *j_it);
	// printf("nb pair %d - %d\n", i, *j_it);
	
	// shortrange, therefore store in simulation.system()
	perturbed_lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
      }
      
    }
  }

  DEBUG(7, "end of function perturbed nonbonded interaction");  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_interaction_spec, typename  t_perturbation_details>
inline void interaction::Perturbed_Nonbonded_Outerloop<
  t_interaction_spec,  t_perturbation_details>
::perturbed_one_four_outerloop(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       Storage & storage)
{
  DEBUG(7, "\tcalculate perturbed 1,4-interactions");
  
  Periodicity_type periodicity(conf.current().box);
  
  std::set<int>::const_iterator it, to;
  std::map<unsigned int, topology::Perturbed_Atom>::const_iterator 
    mit=topo.perturbed_solute().atoms().begin(), 
    mto=topo.perturbed_solute().atoms().end();
  
  for(; mit!=mto; ++mit){
    it = mit->second.one_four_pair().begin();
    to = mit->second.one_four_pair().end();
    
    for( ; it != to; ++it){

      perturbed_one_four_interaction_innerloop
	(topo, conf, mit->second.sequence_number(), *it, periodicity);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_interaction_spec, typename  t_perturbation_details>
inline void interaction::Perturbed_Nonbonded_Outerloop<
  t_interaction_spec,  t_perturbation_details>
::perturbed_RF_excluded_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage)
{

  DEBUG(7, "\tcalculate perturbed excluded RF interactions");

  Periodicity_type periodicity(conf.current().box);

  std::map<unsigned int, topology::Perturbed_Atom>::const_iterator
    mit=topo.perturbed_solute().atoms().begin(),
    mto=topo.perturbed_solute().atoms().end();

  DEBUG(9, "\tSize of perturbed atoms " 
	<< unsigned(topo.perturbed_solute().atoms().size()));
  
  for(; mit!=mto; ++mit){
    perturbed_RF_excluded_interaction_innerloop(topo, conf, mit, periodicity);
  }
}

