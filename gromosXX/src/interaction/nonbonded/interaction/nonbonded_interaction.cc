/**
 * @file nonbonded_interaction.cc
 * template methods of Nonbonded_Interaction.
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
interaction::Nonbonded_Interaction<t_interaction_spec, t_perturbation_spec>
::Nonbonded_Interaction(Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec> *pa)
  : Interaction("NonBonded"),
    Nonbonded_Parameter(),
    m_pairlist_algorithm(pa),
    m_longrange_timing(0.0)
{
}

/**
 * Destructor.
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline 
interaction::Nonbonded_Interaction<t_interaction_spec, t_perturbation_spec>
::~Nonbonded_Interaction()
{
  DEBUG(4, "Nonbonded_Interaction::destructor");
  delete m_pairlist_algorithm;
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline int 
interaction::Nonbonded_Interaction<t_interaction_spec, t_perturbation_spec>
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  const double nonbonded_start = util::now();

  // shared memory do this only once
  m_pairlist_algorithm->prepare(topo, conf, sim);

  typename
    std::vector<Nonbonded_Set<t_interaction_spec, t_perturbation_spec> >::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();

#ifdef OMP
  int tid;
#pragma omp parallel private(tid)
    {
      tid = omp_get_thread_num();
      // calculate the corresponding interactions
      assert(m_nonbonded_set.size() > tid);
      m_nonbonded_set[tid].calculate_interactions(topo, conf, sim,
						  tid, m_omp_num_threads);
    }

#else

  // have to do all from here (probably it's only one, coud unite this,
  // but then maybe it's clearer like it is...
  
  for( ; it != to; ++it){
    it->calculate_interactions(topo, conf, sim);
  }

#endif

  // add the forces, energies, virial...
  it = m_nonbonded_set.begin();

  const double ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & e = conf.current().energies;
  
  for( ; it != to; ++it){
    conf.current().force += it->shortrange_storage().force;

    for(unsigned int i = 0; i < ljs; ++i){
      for(unsigned int j = 0; j < ljs; ++j){
      
	e.lj_energy[i][j] += 
	  it->shortrange_storage().energies.lj_energy[i][j];
	e.crf_energy[i][j] += 
	  it->shortrange_storage().energies.crf_energy[i][j];
      }
    }
    
    if (t_interaction_spec::do_virial){
      DEBUG(7, "\tadd long range virial");

      for(unsigned int i=0; i<3; ++i){
	for(unsigned int j=0; j<3; ++j){

	  DEBUG(8, "set virial = " << it->shortrange_storage().virial_tensor(i,j)
		<< "\tvirial = " << conf.current().virial_tensor(i,j));
	  
	  conf.current().virial_tensor(i,j) +=
	    it->shortrange_storage().virial_tensor(i,j);
	}
      }
      
    }
  }
  
  if (t_perturbation_spec::do_perturbation){

    it = m_nonbonded_set.begin();
      
    for( ; it != to; ++it){
	
      configuration::Energy & pe = conf.current().perturbed_energy_derivatives;
	
      for(unsigned int i = 0; i < ljs; ++i){
	for(unsigned int j = 0; j < ljs; ++j){
      
	  pe.lj_energy[i][j] += 
	    it->shortrange_storage().perturbed_energy_derivatives.lj_energy[i][j];
	  pe.crf_energy[i][j] += 
	    it->shortrange_storage().perturbed_energy_derivatives.crf_energy[i][j];
	}
      }
    } // sets
  } // perturbation

  m_timing += util::now() - nonbonded_start;

  return 0;
  
}

/**
 * calculate the hessian for a given atom.
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
int interaction::Nonbonded_Interaction<t_interaction_spec, t_perturbation_spec>
::calculate_hessian(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    unsigned int atom_i, unsigned int atom_j,
		    math::Matrix & hessian)
{
  typename
    std::vector<Nonbonded_Set<t_interaction_spec, t_perturbation_spec> >::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();

  hessian = 0.0;
  math::Matrix h;

  for( ; it != to; ++it){
    it->calculate_hessian(topo, conf, sim, atom_i, atom_j, h);

    for(unsigned int d1=0; d1 < 3; ++d1)
      for(unsigned int d2=0; d2 < 3; ++d2)
	hessian(d1,d2) += h(d1,d2);
  }
  return 0;
}

/**
 * initialize the arrays
 */
template<typename t_interaction_spec, typename t_perturbation_spec>
inline void interaction::Nonbonded_Interaction<t_interaction_spec, t_perturbation_spec>
::initialize(topology::Topology const & topo,
	     configuration::Configuration const & conf,
	     simulation::Simulation const & sim,
	     bool quiet)
{

#ifdef OMP
  int tid;
#pragma omp parallel private(tid)
    {
      tid = omp_get_thread_num();
      if (tid == 0){
	m_omp_num_threads = omp_get_num_threads();
      }
      
    }
#else
    m_omp_num_threads = 1;
#endif

  DEBUG(15, "nonbonded_interaction::initialize");
  m_nonbonded_set.
    resize(m_omp_num_threads, Nonbonded_Set<t_interaction_spec, t_perturbation_spec>
	      (*this));
  
  typename
    std::vector<Nonbonded_Set<t_interaction_spec, t_perturbation_spec> >::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();
  
  for( ; it != to; ++it){
    it->initialize(topo, conf, sim, quiet);
  }
}


//***************************************************************************
// helper functions 
//***************************************************************************

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void interaction::Nonbonded_Interaction<t_interaction_spec, t_perturbation_spec>
::print_timing(std::ostream & os)
{
      os << "        "
	 << std::setw(36) << std::left << name
	 << std::setw(20) << m_timing << "\n"
	 << "            "
	 << std::setw(32) << std::left << "shortrange"
	 << std::setw(20) 
	 << m_timing - m_pairlist_algorithm->timing()
	 << "\n"
	 << "            "
	 << std::setw(32) << std::left << "longrange"
	  // << std::setw(20) << m_longrange_timing << "\n"
	 << std::setw(20) << "not measured: too expensive" << "\n"

	 << "            "
	 << std::setw(32) << std::left << "pairlist"
	 << std::setw(20) 
	 << m_pairlist_algorithm->timing() - m_longrange_timing<< "\n"
	 << "\n";
      
}
