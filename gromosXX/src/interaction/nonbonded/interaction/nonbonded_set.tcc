/**
 * @file nonbonded_set.tcc
 * template methods of Nonbonded_Set.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline
interaction::Nonbonded_Set<t_interaction_spec, t_perturbation_spec>
::Nonbonded_Set(Nonbonded_Interaction<t_interaction_spec, t_perturbation_spec> & nbi)
  : Nonbonded_Outerloop<t_interaction_spec>(nbi),
    Perturbed_Nonbonded_Outerloop<t_interaction_spec, 
				  typename t_perturbation_spec::perturbation_details>(nbi),
    Perturbed_Nonbonded_Pair<
  t_interaction_spec, typename t_perturbation_spec::perturbation_details>(nbi, 
		      dynamic_cast<Nonbonded_Term&>(*this),
		      dynamic_cast<Perturbed_Nonbonded_Term&>(*this)),
    m_nonbonded_interaction(&nbi)
{
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline int
interaction::Nonbonded_Set<t_interaction_spec, t_perturbation_spec>
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "(Perturbed) Nonbonded_Interaction::calculate_interactions");

  // allow for slow growth (do it every step...)
  if(t_perturbation_spec::do_perturbation)
    set_lambda(topo.lambda(), topo.lambda_exp());

  // zero forces, energies, virial...
  m_shortrange_storage.zero();

  // need to update pairlist?
  if(!(sim.steps() % sim.param().pairlist.skip_step)){
    // create a pairlist

    // double pairlist_start = now();
    
    // zero the longrange forces, energies, virial
    m_longrange_storage.zero();

    DEBUG(7, "\tupdate the parlist");
    m_nonbonded_interaction->pairlist_algorithm().
      update(topo, conf, sim, *this, 0, topo.num_atoms(), 1);
    DEBUG(7, "\tpairlist updated");
    
    // timing.pairlist += now() - pairlist_start;
    // ++timing.count_pairlist;

  }

  // calculate forces / energies
  DEBUG(7, "\tshort range interactions");

  //  double shortrange_start = now();

  lj_crf_outerloop(topo, conf, sim,
		   m_pairlist, m_shortrange_storage);

  if (t_perturbation_spec::do_perturbation){
    DEBUG(7, "\tperturbed short range");
    perturbed_lj_crf_outerloop(topo, conf, sim, 
			       m_perturbed_pairlist,
			       m_shortrange_storage);
  }
  
  // add 1,4 - interactions
  DEBUG(7, "\t1,4 - interactions");
  one_four_outerloop(topo, conf, sim, m_shortrange_storage);
  if(t_perturbation_spec::do_perturbation){
    DEBUG(7, "\tperturbed 1,4 - interactions");
    perturbed_one_four_outerloop(topo, conf, sim, m_shortrange_storage);
  }
  
  // possibly do the RF contributions due to excluded atoms
  if(sim.param().longrange.rf_excluded){
    DEBUG(7, "\tRF excluded interactions and self term");
    RF_excluded_outerloop(topo, conf, sim, m_shortrange_storage);
    if(t_perturbation_spec::do_perturbation){
      DEBUG(7, "\tperturbed RF excluded interactions and self term");
      perturbed_RF_excluded_outerloop(topo, conf, sim, m_shortrange_storage);
    }
  }

  if(t_perturbation_spec::do_perturbation){
    DEBUG(7, "\tperturbed pairs");
    perturbed_pair_outerloop(topo, conf, sim, m_shortrange_storage);
  }

  // add long-range force
  DEBUG(7, "\tadd long range forces");

  m_shortrange_storage.force += m_longrange_storage.force;
  
  // and long-range energies
  DEBUG(7, "\tadd long range energies");
  for(size_t i = 0; i < m_shortrange_storage.energies.lj_energy.size(); ++i){
    for(size_t j = 0; j < m_shortrange_storage.energies.lj_energy.size(); ++j){
      m_shortrange_storage.energies.lj_energy[i][j] += 
	m_longrange_storage.energies.lj_energy[i][j];
      m_shortrange_storage.energies.crf_energy[i][j] += 
	m_longrange_storage.energies.crf_energy[i][j];
    }
  }

  // add longrange virial
  if (t_interaction_spec::do_virial){
    DEBUG(7, "\tadd long range virial");
    for(size_t i=0; i<3; ++i)
      for(size_t j=0; j<3; ++j)
	m_shortrange_storage.virial_tensor(i,j) +=
	  m_longrange_storage.virial_tensor(i,j);
  }
  
  
  if (t_perturbation_spec::do_perturbation){
    // and long-range energy lambda-derivatives
    DEBUG(7, "add long-range lambda-derivatives");

    for(size_t i = 0; 
	i < m_shortrange_storage.perturbed_energy_derivatives.lj_energy.size(); 
	++i){
      for(size_t j = 0; 
	  j < m_shortrange_storage.perturbed_energy_derivatives.lj_energy.size();
	  ++j){

	assert(m_shortrange_storage.perturbed_energy_derivatives.
	       lj_energy.size() > i);
	assert(m_shortrange_storage.perturbed_energy_derivatives.
	       lj_energy[i].size() > j);
	assert(m_shortrange_storage.perturbed_energy_derivatives.
	       lj_energy.size() > j);
	assert(m_shortrange_storage.perturbed_energy_derivatives.
	       lj_energy[j].size() > i);
	
	m_shortrange_storage.perturbed_energy_derivatives.lj_energy[i][j] += 
	  m_longrange_storage.perturbed_energy_derivatives
	  .lj_energy[i][j];
	m_shortrange_storage.perturbed_energy_derivatives.crf_energy[i][j] += 
	  m_longrange_storage.perturbed_energy_derivatives
	  .crf_energy[i][j];
      }
    }
  } // do perturbed

  return 0;
}

/**
 * add a shortrange interaction
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Nonbonded_Set<t_interaction_spec, t_perturbation_spec>
::add_shortrange_pair(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      size_t const i, size_t const j,
		      int pc)
{
  assert(!t_interaction_spec::do_bekker || (pc >= 0 && pc < 27));
  assert(pairlist().size() > i);

  if (t_perturbation_spec::do_perturbation && topo.is_perturbed(i)){

    if (t_perturbation_spec::perturbation_details::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  if(t_interaction_spec::do_bekker)
	    perturbed_pairlist()[i].push_back((pc << 26) + j);
	  else
	    perturbed_pairlist()[i].push_back(j);
	else
	  if(t_interaction_spec::do_bekker)
	    pairlist()[i].push_back((pc << 26) + j);
	  else
	    pairlist()[i].push_back(j);
	return;
      }
    }

    if(t_interaction_spec::do_bekker)
      perturbed_pairlist()[i].push_back((pc << 26) + j);
    else 
      perturbed_pairlist()[i].push_back(j);
  }
  else if (t_perturbation_spec::do_perturbation && topo.is_perturbed(j)){

    if (t_perturbation_spec::perturbation_details::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  if(t_interaction_spec::do_bekker)
	    perturbed_pairlist()[j].push_back(((26 - pc) << 26) + i);
	  else
	    perturbed_pairlist()[j].push_back(i);
	else
	  if(t_interaction_spec::do_bekker)
	    pairlist()[i].push_back((pc << 26) + j);
	  else
	    pairlist()[i].push_back(j);

	return;      
      }
    }
    if(t_interaction_spec::do_bekker)
      perturbed_pairlist()[j].push_back(((26 - pc) << 26) + i);
    else
      perturbed_pairlist()[j].push_back(i);
  }
  else
    if(t_interaction_spec::do_bekker)
      pairlist()[i].push_back((pc << 26) + j);
    else
      pairlist()[i].push_back(j);
     
}

/**
 * add a longrange interaction
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Nonbonded_Set<t_interaction_spec, t_perturbation_spec>
::add_longrange_pair(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     size_t const i, size_t const j,
		     Periodicity_type const & periodicity, int pc)
{
  const double longrange_start = util::now();

  if (t_perturbation_spec::do_perturbation && topo.is_perturbed(i)){

    if (t_perturbation_spec::perturbation_details::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_lj_crf_innerloop(topo, conf, i, j, m_longrange_storage, periodicity, pc);
	else
	  lj_crf_innerloop(topo, conf, i, j, m_longrange_storage, periodicity, pc);
	m_nonbonded_interaction->longrange_timing() += 
	  util::now() - longrange_start;

	return;      
      }
    }
    
    perturbed_lj_crf_innerloop(topo, conf, i, j, m_longrange_storage, periodicity, pc);

  }
  else if (t_perturbation_spec::do_perturbation && topo.is_perturbed(j)){

    if (t_perturbation_spec::perturbation_details::do_scaling){
      // check whether we need to do scaling
      // based on energy groups
      if (sim.param().perturbation.scaled_only){
	std::pair<int, int> 
	  energy_group_pair(topo.atom_energy_group(i),
			    topo.atom_energy_group(j));
	
	if (topo.energy_group_scaling().count(energy_group_pair))
	  perturbed_lj_crf_innerloop(topo, conf, j, i, m_longrange_storage, periodicity, pc);
	else
	  lj_crf_innerloop(topo, conf, i, j, m_longrange_storage, periodicity, pc);

	m_nonbonded_interaction->longrange_timing() += 
	  util::now() - longrange_start;

	return;      
      }
    }
    
    perturbed_lj_crf_innerloop(topo, conf, j, i, m_longrange_storage, periodicity, pc);

  }
  else
    lj_crf_innerloop(topo, conf, i, j, m_longrange_storage, periodicity, pc);

  m_nonbonded_interaction->longrange_timing() += 
    util::now() - longrange_start;

}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Nonbonded_Set<t_interaction_spec, t_perturbation_spec>
::initialize(topology::Topology const & topo,
	     configuration::Configuration const & conf,
	     simulation::Simulation const & sim)
{
  Nonbonded_Outerloop<t_interaction_spec>::initialize(sim);
  if(t_perturbation_spec::do_perturbation)
    Perturbed_Nonbonded_Outerloop<t_interaction_spec, typename t_perturbation_spec::perturbation_details>::initialize(sim);

  m_shortrange_storage.force.resize(conf.current().force.size());
  m_longrange_storage.force.resize(conf.current().force.size());

  m_shortrange_storage.energies.
    resize(conf.current().energies.bond_energy.size(),
	   conf.current().energies.kinetic_energy.size());
  m_longrange_storage.energies.
    resize(conf.current().energies.bond_energy.size(),
	   conf.current().energies.kinetic_energy.size());

  m_shortrange_storage.perturbed_energy_derivatives.resize
    (conf.current().perturbed_energy_derivatives.bond_energy.size(),
     conf.current().perturbed_energy_derivatives.kinetic_energy.size());
  m_longrange_storage.perturbed_energy_derivatives.resize
    (conf.current().perturbed_energy_derivatives.bond_energy.size(),
     conf.current().perturbed_energy_derivatives.kinetic_energy.size());
}

