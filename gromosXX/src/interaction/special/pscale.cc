/**
 * @file pscale.cc
 * contains the implementation
 * for periodic scaling
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * Constructor
 */
template<typename t_interaction_spec>
interaction::Periodic_Scaling<t_interaction_spec>
::Periodic_Scaling(interaction::Forcefield & ff)
  : interaction::Interaction("Periodic_Scaling"),
    m_DI(NULL)
{
  
  // try to get the DihedralAngle interaction from the forcefield
  std::vector<interaction::Interaction *>::iterator
    i_it = ff.begin(),
    i_to = ff.end();

  for( ; i_it != i_to; ++i_it){

    if ((m_DI = dynamic_cast<interaction::Dihedral_Interaction<t_interaction_spec> *>(*i_it)))
      break;
  }
  if (i_it == i_to){
    io::messages.add("Could not access a dihedral interaction.",
		     "Periodic Scaling",
		     io::message::error);
  }

}


/**
 * periodic scaling.
 */
template<typename t_interaction_spec>
int interaction::Periodic_Scaling<t_interaction_spec>
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim)
{
  return 0;
}

/**
 * initialisation
 */
template<typename t_interaction_spec>
int interaction::Periodic_Scaling<t_interaction_spec>
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim,
       bool quiet)
{
  
  if (!quiet)
    std::cout << "PERIODIC SCALING\n";
  
  const int num_J = conf.special().jvalue_curr.size();

  conf.special().pscale.JtoDihedral.resize(num_J, -1);
  conf.special().pscale.KDIH.resize(num_J, 0.0);
  conf.special().pscale.KJ.resize(num_J, 0.0);
  conf.special().pscale.t.resize(num_J, 0.0);
  

  // loop over the J-Values
  std::vector<topology::jvalue_restraint_struct>::const_iterator 
    it = topo.jvalue_restraints().begin(),
    to = topo.jvalue_restraints().end();

    std::cout << std::setw(10) << "dihedral"
	      << std::setw(10) << "i"
	      << std::setw(10) << "j"
	      << std::setw(10) << "k"
	      << std::setw(10) << "l"
	      << std::setw(20) << "type"
	      << std::setw(5) << " "
	      << std::setw(10) << "new type"
	      << "\n";    

  for(int n=0; it != to; ++it, ++n){

    // store force constant (to avoid drift)
    conf.special().pscale.KJ[n] = it->K;

    // look for a matching dihedral angle
    std::vector<topology::four_body_term_struct>::iterator
      d_it = topo.solute().dihedrals().begin(),
      d_to = topo.solute().dihedrals().end();

    for(int d=0; d_it != d_to; ++d_it, ++d){
    
      if (d_it->i == it->i && d_it->j == it->j &&
	  d_it->k == it->k && d_it->l == it->l){
	
	std::cout << std::setw(20) << it->i + 1
		  << std::setw(10) << it->j + 1
		  << std::setw(10) << it->k + 1
		  << std::setw(10) << it->l + 1
		  << std::setw(20) << d_it->type + 1
		  << std::setw(5) << "->";

	// change the type to a new (private) type
	interaction::dihedral_type_struct dts(m_DI->parameter()[d_it->type]);
	m_DI->parameter().push_back(dts);
	d_it->type = m_DI->parameter().size() - 1;

	std::cout << std::setw(10) << d_it->type + 1
		  << "\n";

	// store force constant to avoid drift
	conf.special().pscale.KDIH[n] = m_DI->parameter()[d_it->type].K;
	conf.special().pscale.JtoDihedral[n] = d;
	
      }
    }
    
    
  }

  if (!quiet)
    std::cout << "END\n";

  return 0;
}
