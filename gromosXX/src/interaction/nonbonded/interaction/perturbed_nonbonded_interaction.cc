/**
 * @file perturbed_nonbonded_interaction.cc
 * template methods of Perturbed_Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE nonbonded

#include <util/debug.h>

/**
 * Constructor.
 */
template<typename t_interaction_spec>
inline
interaction::Perturbed_Nonbonded_Interaction<t_interaction_spec>
::Perturbed_Nonbonded_Interaction()
  : Nonbonded_Interaction<t_interaction_spec>(),
    t_interaction_spec::perturbation_filter_type(),
    t_interaction_spec::perturbed_nonbonded_innerloop_type
      (*dynamic_cast<Nonbonded_Base *>(this))
{
}

/**
 * Destructor.
 */
template<typename t_interaction_spec>
inline 
interaction::Perturbed_Nonbonded_Interaction<t_interaction_spec>
::~Perturbed_Nonbonded_Interaction()
{
  DEBUG(4, "Perturbed_Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_interaction_spec>
inline int
interaction::Perturbed_Nonbonded_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "(Perturbed) Nonbonded_Interaction::calculate_interactions");

  // initialize the constants
  if (!sim.steps())
    initialize(topo, conf, sim);

  // allow for slow growth (do it every step...)
  set_lambda(topo.lambda(), topo.lambda_exp());

  // need to update pairlist?
  if(!(sim.steps() % sim.param().pairlist.skip_step)){
    // create a pairlist

    // double pairlist_start = now();
    
    // zero the longrange forces, energies, virial
    force = 0.0;
    energies.zero();
    perturbed_energy_derivatives.zero();
    DEBUG(11,"\tenergy derivative" << perturbed_energy_derivatives.lj_energy[0][0] );
    
    virial_tensor = 0.0;

    DEBUG(7, "\tupdate the parlist");
    m_pairlist_algorithm.update(topo, conf, sim, *this);
    DEBUG(7, "\tpairlist updated");
    
    // timing.pairlist += now() - pairlist_start;
    // ++timing.count_pairlist;

  }

  // calculate forces / energies
  DEBUG(7, "\tshort range interactions");

  //  double shortrange_start = now();

  do_interactions(topo, conf, sim,
		  m_pairlist);
  /*
		  m_pairlist.begin(),
		  m_pairlist.end() );
  */

  if (t_interaction_spec::do_perturbation){
    DEBUG(7, "\tperturbed short range");
    do_perturbed_interactions(topo, conf, sim, 
			      m_perturbed_pairlist.begin(),
			      m_perturbed_pairlist.end() );
  }
  
  // add long-range force
  DEBUG(7, "\tadd long range forces");

  conf.current().force += force;
  
  // and long-range energies
  DEBUG(7, "\tadd long range energies");
  for(size_t i = 0; i < energies.lj_energy.size(); ++i){
    for(size_t j = 0; j < energies.lj_energy.size(); ++j){
      conf.current().energies.lj_energy[i][j] += 
	energies.lj_energy[i][j];
      conf.current().energies.crf_energy[i][j] += 
	energies.crf_energy[i][j];
    }
  }

  // add longrange virial
  if (t_interaction_spec::do_virial){
    DEBUG(7, "\tadd long range virial");
    for(size_t i=0; i<3; ++i)
      for(size_t j=0; j<3; ++j)
	conf.current().virial_tensor(i,j) =
	  conf.current().virial_tensor(i,j) + virial_tensor(i,j);
  }
  
  // add 1,4 - interactions
  DEBUG(7, "\t1,4 - interactions");
  do_14_interactions(topo, conf, sim);
  if(t_interaction_spec::do_perturbation){
    DEBUG(7, "\tperturbed 1,4 - interactions");
    do_perturbed_14_interactions(topo, conf, sim);
  }
  
  // possibly do the RF contributions due to excluded atoms
  if(sim.param().longrange.rf_excluded){
    DEBUG(7, "\tRF excluded interactions and self term");
    do_RF_excluded_interactions(topo, conf, sim);
    if(t_interaction_spec::do_perturbation){
      DEBUG(7, "\tperturbed RF excluded interactions and self term");
      do_perturbed_RF_excluded_interactions(topo, conf, sim);
    }
  }

  if(t_interaction_spec::do_perturbation){
    DEBUG(7, "\tperturbed pairs");
    do_perturbed_pair_interactions(topo, conf, sim);
  }
  
  if (t_interaction_spec::do_perturbation){
    // and long-range energy lambda-derivatives
    DEBUG(7, "add long-range lambda-derivatives");

    for(size_t i = 0; 
	i < perturbed_energy_derivatives.lj_energy.size(); ++i){
      for(size_t j = 0; j < perturbed_energy_derivatives.lj_energy.size(); ++j){

	assert(conf.current().perturbed_energy_derivatives.
	       lj_energy.size() > i);
	assert(conf.current().perturbed_energy_derivatives.
	       lj_energy[i].size() > j);
	assert(conf.current().perturbed_energy_derivatives.
	       lj_energy.size() > j);
	assert(conf.current().perturbed_energy_derivatives.
	       lj_energy[j].size() > i);
	
	conf.current().perturbed_energy_derivatives.lj_energy[i][j] += 
	  perturbed_energy_derivatives
	  .lj_energy[i][j];
	conf.current().perturbed_energy_derivatives.crf_energy[i][j] += 
	  perturbed_energy_derivatives
	  .crf_energy[i][j];
      }
    }
  } // do perturbed

  // timing.shortrange += now() - shortrange_start;
  // ++timing.count_shortrange;

  return 0;
}

/**
 * add a shortrange interaction
 */
template<typename t_interaction_spec>
inline void
interaction::Perturbed_Nonbonded_Interaction<t_interaction_spec>
::add_shortrange_pair(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      size_t const i, size_t const j)
{
  assert(pairlist().size() > i);

  if (perturbed_atom(topo, conf, sim, i)){

    if (t_interaction_spec::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
      
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_pairlist()[i].push_back(j);
	else
	  pairlist()[i].push_back(j);
	return;      
      }
    }
    perturbed_pairlist()[i].push_back(j);
  }
  else if (perturbed_atom(topo, conf, sim, j)){

    if (t_interaction_spec::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_pairlist()[j].push_back(i);
	else
	  pairlist()[i].push_back(j);
	return;      
      }
    }
    
    
    perturbed_pairlist()[j].push_back(i);
  }
  else
    pairlist()[i].push_back(j);
}

/**
 * add a shortrange interaction
 */
template<typename t_interaction_spec>
inline void
interaction::Perturbed_Nonbonded_Interaction<t_interaction_spec>
::add_shortrange_pair(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      size_t const i, size_t const j,
		      int pc)
{
  assert(t_interaction_spec::do_bekker);
  assert(pairlist().size() > i);

  if (perturbed_atom(topo, conf, sim, i)){

    if (t_interaction_spec::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_pairlist()[i].push_back((pc << 26) + j);
	else
	  pairlist()[i].push_back((pc << 26) + j);
	return;
      }
    }

    perturbed_pairlist()[i].push_back((pc << 26) + j);
  }
  else if (perturbed_atom(topo, conf, sim, j)){

    if (t_interaction_spec::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_pairlist()[j].push_back(((26 - pc) << 26) + i);
	else
	  pairlist()[i].push_back((pc << 26) + j);
	return;      
      }
    }
    
    perturbed_pairlist()[j].push_back(((26 - pc) << 26) + i);
  }
  else
    pairlist()[i].push_back((pc << 26) + j);
}

/**
 * add a longrange interaction
 */
template<typename t_interaction_spec>
inline void
interaction::Perturbed_Nonbonded_Interaction<t_interaction_spec>
::add_longrange_pair(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     size_t const i, size_t const j,
		     Periodicity_type const & periodicity)
{
  if (perturbed_atom(topo, conf, sim, i)){

    if (t_interaction_spec::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_interaction_innerloop(topo, conf, i, j, *this, periodicity);
	else
	  interaction_innerloop(topo, conf, i, j, *this, periodicity);
	return;      
      }
    }
    
    perturbed_interaction_innerloop(topo, conf, i, j, *this, periodicity);

  }
  else if (perturbed_atom(topo, conf, sim, j)){

    if (t_interaction_spec::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_interaction_innerloop(topo, conf, j, i, *this, periodicity);
	else
	  interaction_innerloop(topo, conf, i, j, *this, periodicity);
	return;      
      }
    }
    
    perturbed_interaction_innerloop(topo, conf, j, i, *this, periodicity);

  }
  else
    interaction_innerloop(topo, conf, i, j, *this, periodicity);
}

/**
 * add a longrange interaction
 */
template<typename t_interaction_spec>
inline void
interaction::Perturbed_Nonbonded_Interaction<t_interaction_spec>
::add_longrange_pair(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     size_t const i, size_t const j,
		     Periodicity_type const & periodicity, int pc)
{
  if (perturbed_atom(topo, conf, sim, i)){

    if (t_interaction_spec::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_interaction_innerloop(topo, conf, i, j, *this, periodicity, pc);
	else
	  interaction_innerloop(topo, conf, i, j, *this, periodicity, pc);
	return;      
      }
    }
    
    perturbed_interaction_innerloop(topo, conf, i, j, *this, periodicity, pc);

  }
  else if (perturbed_atom(topo, conf, sim, j)){

    if (t_interaction_spec::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_interaction_innerloop(topo, conf, j, i, *this, periodicity, pc);
	else
	  interaction_innerloop(topo, conf, i, j, *this, periodicity, pc);
	return;      
      }
    }
    
    perturbed_interaction_innerloop(topo, conf, j, i, *this, periodicity, pc);

  }
  else
    interaction_innerloop(topo, conf, i, j, *this, periodicity, pc);

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
template<typename t_interaction_spec>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_interaction_spec>
::do_perturbed_interactions(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    Pairlist::iterator it,
			    Pairlist::iterator to)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  Periodicity_type periodicity(conf.current().box);

  if (t_interaction_spec::do_bekker){

    periodicity.recalc_shift_vectors();
    
    int pc;
    size_t j;

    // translate the atom j
    for( ; it != to; ++it){
      
      pc = (*it >> 26);
      j = (*it & 67108863);
      
      perturbed_interaction_innerloop(topo, conf, 
				      it.i(), j, 
				      conf.current(), 
				      periodicity, pc);
    }

  }
  else{
    
    for( ; it != to; ++it){    
      DEBUG(8, "perturbed pair: " << it.i() << " - " << *it);
      perturbed_interaction_innerloop(topo, conf, it.i(), *it, 
				      conf.current(), periodicity);
    }
  }

  DEBUG(7, "end of function perturbed nonbonded interaction");  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_interaction_spec>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_interaction_spec>
::do_perturbed_14_interactions(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim)
{
  DEBUG(7, "\tcalculate perturbed 1,4-interactions");

  Periodicity_type periodicity(conf.current().box);

  std::set<int>::const_iterator it, to;
  std::map<size_t, topology::Perturbed_Atom>::const_iterator 
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
template<typename t_interaction_spec>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_interaction_spec>
::do_perturbed_RF_excluded_interactions(topology::Topology & topo,
					configuration::Configuration & conf,
					simulation::Simulation & sim)
{

  DEBUG(7, "\tcalculate perturbed excluded RF interactions");

  Periodicity_type periodicity(conf.current().box);

  std::map<size_t, topology::Perturbed_Atom>::const_iterator
    mit=topo.perturbed_solute().atoms().begin(),
    mto=topo.perturbed_solute().atoms().end();

  DEBUG(9, "\tSize of perturbed atoms " 
	<< topo.perturbed_solute().atoms().size());
  
  for(; mit!=mto; ++mit){
    perturbed_RF_excluded_interaction_innerloop(topo, conf, mit, periodicity);
  }
}

/**
 * calculate the interactions for the
 * PERTURBED PAIRS
 * (different interaction types in A and in B)
 */
template<typename t_interaction_spec>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_interaction_spec>
::do_perturbed_pair_interactions(topology::Topology & topo,
				 configuration::Configuration & conf,
				 simulation::Simulation & sim)
{
  DEBUG(8, "perturbed pairs");
  
  Periodicity_type periodicity(conf.current().box);
  
  std::vector<topology::perturbed_two_body_term_struct>::const_iterator
    it = topo.perturbed_solute().atompairs().begin(),
    to = topo.perturbed_solute().atompairs().end();
    
  for(; it != to; ++it){
    perturbed_pair_interaction_innerloop(topo, conf, sim, it, periodicity);
  }
  
}

