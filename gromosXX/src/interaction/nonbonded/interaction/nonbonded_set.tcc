/**
 * @file nonbonded_set.tcc
 * template methods of Nonbonded_Set.
 */

// just testing
// the sleep function
#include <unistd.h>

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
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
int
interaction::Nonbonded_Set<t_interaction_spec, t_perturbation_spec>
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 int tid, int num_threads)
{
  DEBUG(4, "(Perturbed) Nonbonded_Interaction::calculate_interactions");

  // allow for slow growth (do it every step...)
  if(t_perturbation_spec::do_perturbation)
    set_lambda(topo.lambda(), topo.lambda_exp());

  // zero forces, energies, virial...
  m_shortrange_storage.zero();

  // need to update pairlist?
  if(!(sim.steps() % sim.param().pairlist.skip_step)){
    //====================
    // create a pairlist
    //====================
    
    // zero the longrange forces, energies, virial
    m_longrange_storage.zero();

    // int atoms_per_thread = topo.num_atoms() / num_threads;
    // int start_atom = atoms_per_thread * tid;
    // int end_atom = atoms_per_thread * (tid + 1);
    
    // m_nonbonded_interaction->pairlist_algorithm().
    // update(topo, conf, sim, *this, start_atom, end_atom, 1);

    // other option: (using the stride)
    // chargegroup based pairlist can only use this one!!!!
    // TODO:
    // move decision to pairlist!!!
    m_nonbonded_interaction->pairlist_algorithm().
      update(topo, conf, sim, *this, tid, topo.num_atoms(), num_threads);

    /*
    sleep(2*tid);
    
    std::cout << "PRINTING OUT THE PAIRLIST\n\n";
    for(size_t i=0; i<100; ++i){
      if (i >= pairlist().size()) break;

      std::cout << "\n\n--------------------------------------------------";
      std::cout << "\n" << i;
      for(size_t j=0; j<pairlist()[i].size(); ++j){

	if (j % 10 == 0) std::cout << "\n\t";
	std::cout << std::setw(7) << pairlist()[i][j];
      }
    }
    */
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
  if (tid == 0){
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

    // one we definitely have...
    const size_t lj_size 
      = m_shortrange_storage.perturbed_energy_derivatives[0].lj_energy.size();
    
    DEBUG(8, "from a set of "
	  << m_shortrange_storage.perturbed_energy_derivatives.size()
	  << " lambda deps");
    
    for(size_t s=0, s_to = m_shortrange_storage.perturbed_energy_derivatives.size();
	s != s_to; ++s){
      for(size_t i = 0; i < lj_size; ++i){
	for(size_t j = 0; j < lj_size; ++j){

	  assert(m_shortrange_storage.perturbed_energy_derivatives[s].
		 lj_energy.size() > i);
	  assert(m_shortrange_storage.perturbed_energy_derivatives[s].
		 lj_energy[i].size() > j);
	  assert(m_shortrange_storage.perturbed_energy_derivatives[s].
		 lj_energy.size() > j);
	  assert(m_shortrange_storage.perturbed_energy_derivatives[s].
		 lj_energy[j].size() > i);
	
	  m_shortrange_storage.perturbed_energy_derivatives[s].lj_energy[i][j] += 
	    m_longrange_storage.perturbed_energy_derivatives[s]
	    .lj_energy[i][j];
	  m_shortrange_storage.perturbed_energy_derivatives[s].crf_energy[i][j] += 
	    m_longrange_storage.perturbed_energy_derivatives[s]
	    .crf_energy[i][j];
	}
      }
    }
  } // do perturbed

  return 0;
}

/**
 * calculate the hessian for a given atom.
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline int
interaction::Nonbonded_Set<t_interaction_spec, t_perturbation_spec>
::calculate_hessian(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    size_t const atom_i, size_t const atom_j,
		    math::Matrix & hessian){
  
  hessian = 0.0;
  
  // loop over the pairlist

  //*************************
  // standard implementation
  //*************************

  std::vector<size_t>::const_iterator j_it, j_to;
  size_t i;
  size_t size_i = pairlist().size();
  math::Vec r;
  math::Matrix h;
  Periodicity_type periodicity(conf.current().box);
  
  for(j_it = pairlist()[atom_i].begin(),
	j_to = pairlist()[i].end();
      j_it != j_to;
      ++j_it){

    if (*j_it == atom_j){
      periodicity.nearest_image(conf.current().pos(atom_i),
				conf.current().pos(atom_j),
				r);
    }
    else continue;
      
    const lj_parameter_struct &lj = 
      m_nonbonded_interaction->lj_parameter(topo.iac(atom_i),
					    topo.iac(atom_j));
    
    lj_crf_hessian(r,
		   lj.c6, lj.c12, 
		   topo.charge()(atom_i) * topo.charge()(atom_j),
		   h);

    for(size_t d1=0; d1 < 3; ++d1)
      for(size_t d2=0; d2 < 3; ++d2)
	hessian(d1,d2) += h(d1,d2);
  }
  // and the other way round
  for(j_it = pairlist()[atom_j].begin(),
	j_to = pairlist()[i].end();
      j_it != j_to;
      ++j_it){

    if (*j_it == atom_i){
      periodicity.nearest_image(conf.current().pos(atom_i),
				conf.current().pos(atom_j),
				r);
    }
    else continue;
      
    const lj_parameter_struct &lj = 
      m_nonbonded_interaction->lj_parameter(topo.iac(atom_i),
					    topo.iac(atom_j));
    
    lj_crf_hessian(r,
		   lj.c6, lj.c12, 
		   topo.charge()(atom_i) * topo.charge()(atom_j),
		   h);

    for(size_t d1=0; d1 < 3; ++d1)
      for(size_t d2=0; d2 < 3; ++d2)
	hessian(d1,d2) += h(d1,d2);
  }

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
  DEBUG(8, "add shortrange pair " << i << " - " << j);

  assert(!t_interaction_spec::do_bekker || (pc >= 0 && pc < 27));
  assert(pairlist().size() > i);

  if (t_perturbation_spec::do_perturbation){
    if (topo.is_perturbed(i)){
      
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
    else if (topo.is_perturbed(j)){

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
    else{
      if(t_interaction_spec::do_bekker)
	pairlist()[i].push_back((pc << 26) + j);
      else
	pairlist()[i].push_back(j);
    }
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
  DEBUG(8, "add longrange pair " << i << " - " << j);

  // const double longrange_start = util::now();

  if (t_perturbation_spec::do_perturbation){
    if (topo.is_perturbed(i)){

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
	  // m_nonbonded_interaction->longrange_timing() += 
	  // util::now() - longrange_start;
	  
	  return;      
	}
      }
      
      perturbed_lj_crf_innerloop(topo, conf, i, j, m_longrange_storage, periodicity, pc);
      
    }
    else if (topo.is_perturbed(j)){
      
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
	  
	  // m_nonbonded_interaction->longrange_timing() += 
	  // util::now() - longrange_start;
	  
	  return;      
	}
      }
      
      perturbed_lj_crf_innerloop(topo, conf, j, i, m_longrange_storage, periodicity, pc);
      
    }
    else
      lj_crf_innerloop(topo, conf, i, j, m_longrange_storage, periodicity, pc);
  }
  
  else
    lj_crf_innerloop(topo, conf, i, j, m_longrange_storage, periodicity, pc);

  // m_nonbonded_interaction->longrange_timing() += 
  // util::now() - longrange_start;
  
}

template<typename t_interaction_spec, typename t_perturbation_spec>
void
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

  
  size_t es = conf.current().perturbed_energy_derivatives.size();
  
  DEBUG(7, "initialize lambda dep to " << es);

  m_shortrange_storage.perturbed_energy_derivatives.
    resize(es);
  m_longrange_storage.perturbed_energy_derivatives.
    resize(es);

  for(size_t s = 0; s < es; ++s){

    m_shortrange_storage.perturbed_energy_derivatives[s].resize
      (conf.current().perturbed_energy_derivatives[s].bond_energy.size(),
       conf.current().perturbed_energy_derivatives[s].kinetic_energy.size());
    m_longrange_storage.perturbed_energy_derivatives[s].resize
      (conf.current().perturbed_energy_derivatives[s].bond_energy.size(),
       conf.current().perturbed_energy_derivatives[s].kinetic_energy.size());
  }
  
  // and the pairlists
  pairlist().resize(topo.num_atoms());
  if(t_perturbation_spec::do_perturbation){
    perturbed_pairlist().resize(topo.num_atoms());
  }

  // check if we can guess the number of pairs
  const double vol = math::volume(conf.current().box, conf.boundary_type);
  if (vol){
    const double c3 = sim.param().pairlist.cutoff_short *
      sim.param().pairlist.cutoff_short *
      sim.param().pairlist.cutoff_short;
    
    const size_t pairs = 
      int(1.3 * topo.num_atoms() / vol * 4.0 / 3.0 * math::Pi * c3);

    std::cout << "\n\testimated pairlist size (per atom) : "
	      << pairs << "\n\n";
    
    for(size_t i=0; i<topo.num_atoms(); ++i)
      pairlist()[i].reserve(pairs);
    
  }

}

