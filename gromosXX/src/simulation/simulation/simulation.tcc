/**
 * @file simulation.tcc
 * contains the (template) methods
 * for the class Simulation
 */

/**
 * Constructor
 */
template<typename t_topo, typename t_system>
inline simulation::Simulation<t_topo, t_system>
::Simulation(t_topo &topo, t_system &sys)
  : m_topology(topo),
    m_system(sys),
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
  system().energies().kinetic_energy.resize(multibath().size());
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
