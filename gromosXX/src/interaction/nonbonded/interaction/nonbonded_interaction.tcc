/**
 * @file nonbonded_interaction.tcc
 * template methods of Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

/**
 * Constructor.
 */
template<typename t_interaction_spec, bool perturbed>
inline
interaction::Nonbonded_Interaction<t_interaction_spec, perturbed>
::Nonbonded_Interaction(Pairlist_Algorithm<t_interaction_spec, perturbed> *pa)
  : Interaction("NonBonded"),
    Nonbonded_Parameter(),
    m_pairlist_algorithm(pa)
{
}

/**
 * Destructor.
 */
template<typename t_interaction_spec, bool perturbed>
inline 
interaction::Nonbonded_Interaction<t_interaction_spec, perturbed>
::~Nonbonded_Interaction()
{
  DEBUG(4, "Nonbonded_Interaction::destructor");
  delete m_pairlist_algorithm;
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_interaction_spec, bool perturbed>
inline int 
interaction::Nonbonded_Interaction<t_interaction_spec, perturbed>
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  const double nonbonded_start = util::now();

  m_pairlist_algorithm->prepare();

  typename
    std::vector<Nonbonded_Set<t_interaction_spec, perturbed> >::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();
  
  for( ; it != to; ++it){
    it->calculate_interactions(topo, conf, sim);
  }

  // add the forces, energies, virial...
  it = m_nonbonded_set.begin();

  const double ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & e = conf.current().energies;
  
  for( ; it != to; ++it){
    conf.current().force += it->shortrange_storage().force;

    for(size_t i = 0; i < ljs; ++i){
      for(size_t j = 0; j < ljs; ++j){
      
	e.lj_energy[i][j] += 
	  it->shortrange_storage().energies.lj_energy[i][j];
	e.crf_energy[i][j] += 
	  it->shortrange_storage().energies.crf_energy[i][j];
      }
    }
    
    if (t_interaction_spec::do_virial){
      DEBUG(7, "\tadd long range virial");

      for(size_t i=0; i<3; ++i)
	for(size_t j=0; j<3; ++j)
	  conf.current().virial_tensor(i,j) +=
	    it->shortrange_storage().virial_tensor(i,j);
    }
  }
  
  if (perturbed){

    it = m_nonbonded_set.begin();

    configuration::Energy & pe = conf.current().perturbed_energy_derivatives;
  
    for( ; it != to; ++it){

      for(size_t i = 0; i < ljs; ++i){
	for(size_t j = 0; j < ljs; ++j){
      
	  pe.lj_energy[i][j] += 
	    it->shortrange_storage().perturbed_energy_derivatives.lj_energy[i][j];
	  pe.crf_energy[i][j] += 
	    it->shortrange_storage().perturbed_energy_derivatives.crf_energy[i][j];
	}
      }
    }
  }

  m_timing += util::now() - nonbonded_start;

  return 0;
  
}

/**
 * initialize the arrays
 */
template<typename t_interaction_spec, bool perturbed>
inline void interaction::Nonbonded_Interaction<t_interaction_spec, perturbed>
::initialize(topology::Topology const & topo,
	     configuration::Configuration const & conf,
	     simulation::Simulation const & sim)
{
  DEBUG(15, "nonbonded_interaction::initialize");
  m_nonbonded_set.
    resize(1,Nonbonded_Set<t_interaction_spec, perturbed>
	      (*this));
  
  typename
    std::vector<Nonbonded_Set<t_interaction_spec, perturbed> >::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();
  
  for( ; it != to; ++it){
    it->initialize(topo, conf, sim);
  }
}


//***************************************************************************
// helper functions 
//***************************************************************************

template<typename t_interaction_spec, bool perturbed>
inline void interaction::Nonbonded_Interaction<t_interaction_spec, perturbed>
::print_timing(std::ostream & os)
{
      os << "        "
	 << std::setw(36) << std::left << name
	 << std::setw(20) << m_timing << "\n"
	 << "            "
	 << std::setw(32) << std::left << "shortrange"
	//	 << std::setw(20) << m_shortrange_timing << "\n"
	 << "            "
	 << std::setw(32) << std::left << "pairlist/longrange"
	// << std::setw(20) << m_pairlist_timing << "\n"
	 << "                "
	 << std::setw(28) << std::left << "longrange forces"
	//	 << std::setw(20) << m_longrange_timing 
	 << "\n";
      
}

