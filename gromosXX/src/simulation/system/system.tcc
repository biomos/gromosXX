/**
 * @file system.tcc
 * inline methods definition
 */

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE system

#include "../../debug.h"

/**
 * Constructor
 */
template<math::boundary_enum b>
inline simulation::System<b>::System()
  : m_position1(0),
    m_position2(0),
    m_velocity1(0),
    m_velocity2(0),
    m_force1(0),
    m_force2(0),
    m_periodicity()
{
  m_pos = &m_position1;
  m_old_pos = &m_position2;
  m_vel = &m_velocity1;
  m_old_vel = &m_velocity2;
  m_force = &m_force1;
  m_old_force = &m_force2;
}

/**
 * set the number of atoms.
 * using resizeAndPreserve. Therefore
 * you can enlarge the system (or shrink it)
 * while keeping all existing positions/velocities/...
 * a faster version would be just resize, but then
 * the arrays contain garbage...
 */
template<math::boundary_enum b>
inline void simulation::System<b>::resize(size_t s)
{
  DEBUG(7, "system resize: " << s);
  
  m_position1.resizeAndPreserve(s);
  m_position2.resizeAndPreserve(s);
  m_velocity1.resizeAndPreserve(s);
  m_velocity2.resizeAndPreserve(s);
  m_force1.resizeAndPreserve(s);
  m_force2.resizeAndPreserve(s);
}

/**
 * position accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::pos()
{
  return *m_pos;
}
/**
 * old position accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::old_pos()const
{
  return *m_old_pos;
}

/**
 * velocity accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::vel()
{
  return *m_vel;
}

/**
 * const old velocity accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::old_vel()const
{
  return *m_old_vel;
}

/**
 * old velocity accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::old_vel()
{
  return *m_old_vel;
}

/**
 * force accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::force()
{
  return *m_force;
}

/**
 * old force accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::old_force()const
{
  return *m_old_force;
}

/**
 * exchange positions
 */
template<math::boundary_enum b>
inline void simulation::System<b>::exchange_pos()
{
  math::VArray *dummy = m_pos;
  m_pos = m_old_pos;
  m_old_pos = dummy;
}
/**
 * exchange velocities
 */
template<math::boundary_enum b>
inline void simulation::System<b>::exchange_vel()
{
  math::VArray *dummy = m_vel;
  m_vel = m_old_vel;
  m_old_vel = dummy;
}
/**
 * exchange forces
 */
template<math::boundary_enum b>
inline void simulation::System<b>::exchange_force()
{
  math::VArray *dummy = m_force;
  m_force = m_old_force;
  m_old_force = dummy;
}

/**
 * periodicity accessor.
 */
template<math::boundary_enum b>
inline math::Periodicity<b> & simulation::System<b>
::periodicity()
{
  return m_periodicity;
}

/**
 * const periodicity accessor.
 */
template<math::boundary_enum b>
inline math::Periodicity<b> const & simulation::System<b>
::periodicity()const
{
  return m_periodicity;
}

/**
 * box index accessor.
 */
template<math::boundary_enum b>
inline std::vector<int[3]> & simulation::System<b>::box_index()
{
  return m_box_index;
}

template<math::boundary_enum b>
inline int[3] & simulation::System<b>::box_index(size_t i)
{
  assert(i < m_box_index.size());
  return m_box_index[i];
}

namespace simulation
{
  template<math::boundary_enum b>
  inline std::ostream &operator<<(std::ostream &os, System<b> &sys)
  {
    os << "a system";
    return os;
  }
}

