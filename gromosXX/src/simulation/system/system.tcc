/**
 * @file system.tcc
 * inline methods definition
 */

/**
 * Constructor
 */
inline simulation::system::system()
  : m_position1(0),
    m_position2(0),
    m_velocity1(0),
    m_velocity2(0),
    m_force1(0),
    m_force2(0)
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
inline void simulation::system::resize(size_t s)
{
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
inline math::VArray & simulation::system::pos()
{
  return *m_pos;
}
/**
 * old position accessor
 */
inline math::VArray const & simulation::system::old_pos()const
{
  return *m_old_pos;
}

/**
 * velocity accessor
 */
inline math::VArray & simulation::system::vel()
{
  return *m_vel;
}

/**
 * old velocity accessor
 */
inline math::VArray const & simulation::system::old_vel()const
{
  return *m_old_vel;
}

/**
 * force accessor
 */
inline math::VArray & simulation::system::force()
{
  return *m_force;
}

/**
 * old force accessor
 */
inline math::VArray const & simulation::system::old_force()const
{
  return *m_old_force;
}

/**
 * exchange positions
 */
inline void simulation::system::exchange_pos()
{
  math::VArray *dummy = m_pos;
  m_pos = m_old_pos;
  m_old_pos = dummy;
}
/**
 * exchange velocities
 */
inline void simulation::system::exchange_vel()
{
  math::VArray *dummy = m_vel;
  m_vel = m_old_vel;
  m_old_vel = dummy;
}
/**
 * exchange forces
 */
inline void simulation::system::exchange_force()
{
  math::VArray *dummy = m_force;
  m_force = m_old_force;
  m_old_force = dummy;
}
/**
 * box accessor.
 */
inline math::Matrix & simulation::system::box()
{
  return m_box;
}

/**
 * boundary condition accessor.
 */
inline simulation::boundary_enum simulation::system::boundary_condition()
{
  return m_boundary_condition;
}

/**
 * set boundary condition.
 */
inline void simulation::system::boundary_condition(boundary_enum b)
{
  m_boundary_condition = b;
}

namespace simulation
{
  inline std::ostream &operator<<(std::ostream &os, system &sys)
  {
    os << "a system";
    return os;
  }
}

