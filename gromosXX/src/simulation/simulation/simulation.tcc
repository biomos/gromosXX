/**
 * @file simulation.tcc
 * contains the (template) methods
 * for the class simulation
 */

/**
 * Constructor
 */
template<typename t_topo, typename t_system>
inline simulation::simulation<t_topo, t_system>
::simulation(t_topo &topo, t_system &sys)
  : m_topology(topo),
    m_system(sys)
{
}

/**
 * topology accessor
 */
template<typename t_topo, typename t_system>
inline t_topo & simulation::simulation<t_topo, t_system>
::topology()
{
  return m_topology;
}

/**
 * system accessor
 */
template<typename t_topo, typename t_system>
inline t_system & simulation::simulation<t_topo, t_system>
::system()
{
  return m_system;
}

namespace simulation
{
  template<typename t_topo, typename t_system>
  std::ostream &operator<<(std::ostream &os, simulation<t_topo, t_system> &sys)
  {
    os << "a simulation";
    return os;
  }
}
