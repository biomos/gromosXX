/**
 * @file simulation.tcc
 * contains the (template) methods
 * for the class Simulation
 */

// BZ_USING_NAMESPACE(blitz)

/**
 * Constructor
 */
template<typename t_topo, typename t_system>
inline simulation::Simulation<t_topo, t_system>
::Simulation()
  : m_topology(),
    m_system(),
    m_time(0.0),
    m_old_time(0.0),
    m_steps(0)
{
}

/**
 * const topology accessor
 */
template<typename t_topo, typename t_system>
inline t_topo const & simulation::Simulation<t_topo, t_system>
::topology()const
{
  return m_topology;
}

/**
 * topology accessor
 */
template<typename t_topo, typename t_system>
inline t_topo & simulation::Simulation<t_topo, t_system>
::topology()
{
  return m_topology;
}

/**
 * const system accessor
 */
template<typename t_topo, typename t_system>
inline t_system const & simulation::Simulation<t_topo, t_system>
::system()const
{
  return m_system;
}

/**
 * system accessor
 */
template<typename t_topo, typename t_system>
inline t_system & simulation::Simulation<t_topo, t_system>
::system()
{
  return m_system;
}

/**
 * time accessor
 */
template<typename t_topo, typename t_system>
inline double simulation::Simulation<t_topo, t_system>
::time()
{
  return m_time;
}

/**
 * old time accessor
 */
template<typename t_topo, typename t_system>
inline double simulation::Simulation<t_topo, t_system>
::old_time()
{
  return m_old_time;
}

/**
 * set time
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>
::time(double t)
{
  m_time = t;
}

/**
 * increase time by dt
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>
::increase_time(double dt)
{
  m_old_time=m_time;
  m_time += dt;
  ++m_steps;
}

/**
 * steps accessor
 */
template<typename t_topo, typename t_system>
inline int simulation::Simulation<t_topo, t_system>
::steps()
{
  return m_steps;
}

/**
 * Nonbonded interaction class
 */
template<typename t_topo, typename t_system>
inline simulation::Nonbonded & simulation::Simulation<t_topo, t_system>
::nonbonded()
{
  return m_nonbonded;
}

/**
 * const nonbonded interaction class
 */
template<typename t_topo, typename t_system>
inline simulation::Nonbonded const & simulation::Simulation<t_topo, t_system>
::nonbonded()const
{
  return m_nonbonded;
}

/**
 * multibath parameter
 */
template<typename t_topo, typename t_system>
inline simulation::Multibath const & simulation::Simulation<t_topo, t_system>
::multibath()const
{
  return m_multibath;
}

/**
 * multibath parameter
 */
template<typename t_topo, typename t_system>
inline simulation::Multibath & simulation::Simulation<t_topo, t_system>
::multibath()
{
  return m_multibath;
}

/**
 * add solvent molecules to the simulation (system).
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>::solvate(size_t solv, size_t num_molecules)
{
  topology().solvate(solv, num_molecules);
  system().resize(topology().num_solute_atoms() + 
		  topology().num_solvent_atoms());
}

/**
 * calculate the degrees of freedom.
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>::calculate_degrees_of_freedom()
{
  multibath().calculate_degrees_of_freedom(topology());
  system().energies().resize(0, multibath().size());
}

/**
 * put chargegroups into box while updating the box indices for the atoms.
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>::put_chargegroups_into_box()
{
  io::messages.add("not implemented","Simulation: put_chargegroups_into_box",
		   io::message::error);
}

/**
 * pressure calculation
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>
::pressure_calculation(bool pc)
{
  m_pressure_calculation = pc;
}

/**
 * pressure calculation
 */
template<typename t_topo, typename t_system>
inline bool simulation::Simulation<t_topo, t_system>
::pressure_calculation()
{
  return m_pressure_calculation;
}
/**
 * pressure calculation
 */
template<typename t_topo, typename t_system>
inline const bool simulation::Simulation<t_topo, t_system>
::pressure_calculation()const
{
  return m_pressure_calculation;
}


/**
 * calculate positions relative to molecular com
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>::
calculate_mol_com()
{  
  Molecule_Iterator m_it = topology().molecule_begin(),
    m_to = topology().molecule_end();
  
  math::Vec com_pos;
  math::Matrix com_ekin;
  
  system().molecular_kinetic_energy() = 0.0;

  for( ; m_it != m_to; ++m_it){
    system().center_of_mass(m_it.begin(), m_it.end(),
				topology().mass(),
				com_pos, com_ekin);

    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	system().molecular_kinetic_energy()(i,j) += com_ekin(i,j);
    
    Atom_Iterator a_it = m_it.begin(),
      a_to = m_it.end();
    
    math::VArray &pos = system().pos();

    for( ; a_it != a_to; ++a_it){
      assert(unsigned(system().rel_mol_com_pos().size()) > *a_it);
      system().periodicity().nearest_image(pos(*a_it), com_pos, 
					   system().rel_mol_com_pos()(*a_it));
    }

  }
}

/**
 * ekin term for atomic virial
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>::
calculate_atom_ekin()
{  
  system().molecular_kinetic_energy() = 0.0;

  for(size_t i=0; i < topology().num_atoms(); ++i){
    for(int a=0; a<3; ++a){
      for(int b=0; b<3; ++b){
	system().molecular_kinetic_energy()(a, b) +=
	  0.5 * topology().mass()(i) * 
	  system().vel()(i)(a) * system().vel()(i)(b);
      }
    }
  }

  // system().molecular_kinetic_energy() *= 0.5;
}

/**
 * calculate molecular kinetic energies
 * (translational and rotational+internal)
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>::
calculate_mol_ekin()
{  
  Molecule_Iterator m_it = topology().molecule_begin(),
    m_to = topology().molecule_end();
  
  math::Vec com_v, new_com_v;
  double com_ekin, ekin, new_com_ekin, new_ekin;
  
  // assume zero energies has been called!
  // but not in the multibath structs...
  for(size_t i=0; i < multibath().size(); ++i)
    multibath().bath(i).ekin = 0.0;

  size_t ir_bath, com_bath;
  
  for( ; m_it != m_to; ++m_it){
    system().molecular_translational_ekin(m_it.begin(), m_it.end(),
					  topology().mass(),
					  com_v, com_ekin, ekin,
					  new_com_v, new_com_ekin, new_ekin);

    
    multibath().in_bath(*m_it.begin(), com_bath, ir_bath);

    DEBUG(9, "adding to bath: com: "
	  << com_bath << " ir: " << ir_bath);
    DEBUG(10, "number of baths: energy " 
	  << system().energies().kinetic_energy.size()
	  << " com ekin "
	  << system().energies().com_kinetic_energy.size()
	  << " ir ekin "
	  << system().energies().ir_kinetic_energy.size()
	  );
    
    // store the new ones in multibath (velocity scaling for next step)
    multibath().bath(com_bath).ekin += new_com_ekin;
    multibath().bath(ir_bath).ekin += new_ekin - new_com_ekin;

    // and the averages in the energies
    system().energies().com_kinetic_energy[com_bath] += com_ekin;
    system().energies().ir_kinetic_energy[ir_bath] += ekin - com_ekin;

  }

  // loop over the bath kinetic energies
  for(size_t i=0; i<system().energies().kinetic_energy.size(); ++i)
    system().energies().kinetic_energy[i] =
      system().energies().com_kinetic_energy[i] +
      system().energies().ir_kinetic_energy[i];

}

template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>::remove_com_motion
  (double const dt, bool remove_trans, bool remove_rot,
   double &ekin_trans, double &ekin_rot)
{
#ifndef NDEBUG
  std::cout.precision(20);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
#endif

  math::Vec com_v = 0.0, com_r = 0.0;
  double com_mass = 0.0;
  
  for(size_t i = 0; i < topology().num_atoms(); ++i){
    com_mass += topology().mass()(i);
    com_v += topology().mass()(i) * system().vel()(i);
    com_r += topology().mass()(i) * system().pos()(i) - 
      0.5 * topology().mass()(i) * system().vel()(i)*dt;
  }
  com_v /= com_mass;
  com_r /= com_mass;

  DEBUG(7, "totmass " << com_mass);
  DEBUG(7, "com_v " << com_v);
  DEBUG(7, "com_r " << com_r);
  DEBUG(7, "com_Ekin " << 0.5*com_mass*dot(com_v,com_v));

  ekin_trans = 0.5*com_mass*dot(com_v,com_v);
  
  math::Vec com_L = 0.0;
  math::Matrix com_I;
  com_I.initialize(0.0);
  
  for(size_t i = 0; i < topology().num_atoms(); ++i){
    math::Vec r = system().pos()(i) - 0.5*dt*system().vel()(i) - com_r;
    DEBUG(7, "pos  " << i << " " << system().pos()(i) - 0.5*dt*system().vel()(i));
    DEBUG(7, "posp " << i << " " << r);
    
    //com_L += topology().mass()(i) * 
    //  math::cross(r, system().vel()(i) - com_v);    
    com_L += topology().mass()(i)* math::cross(system().pos()(i), system().vel()(i));

    // inertia tensor
    double r2 = dot(r,r);
    com_I(0,0) += topology().mass()(i) * (r(1)*r(1)+r(2)*r(2));
    com_I(1,1) += topology().mass()(i) * (r(0)*r(0)+r(2)*r(2));
    com_I(2,2) += topology().mass()(i) * (r(0)*r(0)+r(1)*r(1));
    com_I(1,0) += topology().mass()(i) * (-r(0)*r(1));
    com_I(0,1) += topology().mass()(i) * (-r(0)*r(1));
    com_I(2,0) += topology().mass()(i) * (-r(0)*r(2));
    com_I(0,2) += topology().mass()(i) * (-r(0)*r(2));
    com_I(2,1) += topology().mass()(i) * (-r(1)*r(2));
    com_I(1,2) += topology().mass()(i) * (-r(1)*r(2));
    
  }
  com_L -= com_mass*math::cross(com_r, com_v);
  
  DEBUG(7, "Angular momentum " << com_L);
  
  // invert the inertia tensor
  math::Matrix com_II;
  const double denom = -com_I(2,0)*com_I(2,0)*com_I(1,1)
    + 2 * com_I(0,1) * com_I(0,2) * com_I(1,2)
    - com_I(0, 0) * com_I(1,2) * com_I(1,2)
    - com_I(0,1) * com_I(0,1) * com_I(2,2)
    + com_I(0,0) * com_I(1,1) * com_I(2,2);
  
  com_II(0,0) = (-com_I(1,2)*com_I(1,2) + com_I(1,1) * com_I(2,2));
  com_II(1,0) = com_II(0,1) = (com_I(0,2) * com_I(1,2)
    - com_I(0,1) * com_I(2,2));
  com_II(0,2) = com_II(2,0) = (-com_I(0,2)*com_I(1,1)
    + com_I(0,1)*com_I(1,2));

  com_II(1,1) = (-com_I(0,2)*com_I(0,2) + com_I(0,0) * com_I(2,2));
  com_II(1,2) = com_II(2,1) = (com_I(0,1)*com_I(0,2)
    - com_I(0,0) * com_I(1,2));

  com_II(2,2) = (-com_I(0,1)*com_I(0,1) + com_I(0,0)*com_I(1,1));
  
  math::Matrix inv;
  DEBUG(7, "inertia tensor:\n"<< com_I);
  DEBUG(7, "determinant : " << denom);
  DEBUG(7, "inverted tens :\n" << com_II);

#ifndef NDEBUG
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      inv(i,j)=0.0;
      for(int k=0; k<3; k++){
	inv(i,j)+= com_I(i,k)*com_II(k,j);
      }
    }
  }
#endif 

  DEBUG(7, "one matrix\n" << inv);
  
  // get the angular velocity around the COM
  math::Vec com_O;
  com_O = blitz::product(com_II, com_L) / denom;
  
  DEBUG(7, " angular velocity " << com_O);
  
  DEBUG(7, " com_Ekin_rot " << 0.5*dot(com_O,com_L));
 
  ekin_rot =  0.5*dot(com_O,com_L);
  
  if(remove_rot){
    DEBUG(7, "removing rot");
    for(size_t i=0; i<topology().num_atoms(); ++i){
      math::Vec r = system().pos()(i) -0.5*dt*system().vel()(i) - com_r;
      system().vel()(i) -=math::cross(com_O, r); 
    }
  }
  if(remove_trans){
    DEBUG(7, "removing trans");
    
    // get the corrected velocities
    for(size_t i=0; i<topology().num_atoms(); ++i){
      system().vel()(i) -= com_v;
    }
  }

}

template<typename t_topo, typename t_system>
inline int simulation::Simulation<t_topo, t_system>::check_state()const
{
  DEBUG(7, "checking state of Simulation");

  int result = 0;

  // check associated classes
  result += m_topology.check_state();
  result += m_system.check_state();

  result += m_nonbonded.check_state();
  result += m_multibath.check_state(m_topology.num_atoms());

  if (m_time < 0)
    io::messages.add("Simulation time < 0", "Simulation::check_state",
		     io::message::warning);
  if (m_steps < 0){
    io::messages.add("Simulation steps < 0", "Simulation::check_state",
		     io::message::error);
    ++result;
  }

  return result;
  
}

namespace simulation
{
  template<typename t_topo, typename t_system>
  std::ostream &operator<<(std::ostream &os, Simulation<t_topo, t_system> &sys)
  {
    os << "a simulation";
    return os;
  }
}
