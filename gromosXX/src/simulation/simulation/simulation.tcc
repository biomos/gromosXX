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
    m_steps(0),
    m_nonbonded_update(5),
    m_nonbonded_cutoff_short(0.8),
    m_nonbonded_cutoff_long(1.4)
{
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
 * pairlist update every n steps.
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>
::nonbonded_update(int const update_step)
{
  m_nonbonded_update = update_step;
}
  
/**
 * accessor pairlist update.
 */
template<typename t_topo, typename t_system>
inline int simulation::Simulation<t_topo, t_system>
::nonbonded_update()const
{
  return m_nonbonded_update;
}

/**
 * set short range cutoff.
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>
::nonbonded_cutoff_short(double const cutoff_short)
{
  m_nonbonded_cutoff_short = cutoff_short;
}

/**
 * get short range cutoff.
 */
template<typename t_topo, typename t_system>
inline double simulation::Simulation<t_topo, t_system>
::nonbonded_cutoff_short()const
{
  return m_nonbonded_cutoff_short;
}

/**
 * set long range cutoff.
 */
template<typename t_topo, typename t_system>
inline void simulation::Simulation<t_topo, t_system>
::nonbonded_cutoff_long(double const cutoff_long)
{
  m_nonbonded_cutoff_long = cutoff_long;
}

/**
 * get long range cutoff.
 */
template<typename t_topo, typename t_system>
inline double simulation::Simulation<t_topo, t_system>
::nonbonded_cutoff_long()const
{
  return m_nonbonded_cutoff_long;
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
