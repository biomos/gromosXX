/**
 * @file topology.tcc
 * inline methods definition
 */

/**
 * Constructor
 */
inline simulation::topology::topology()
  : m_mass(0)
{
}

/**
 * mass accessor
 */
math::SArray & simulation::topology::mass()
{
  return m_mass;
}

/**
 * const mass accessor
 */
math::SArray const & simulation::topology::mass()const
{
  return m_mass;
}

/**
 * the number of solute atoms
 */
inline size_t simulation::topology::num_solute_atoms()const
{
  return m_solute_atoms;
}

/**
 * set the number of solute atoms (and resize
 * the apropriate arrays.
 */
inline void simulation::topology::num_solute_atoms(size_t atoms)
{
  m_solute_atoms = atoms;
  m_mass.resizeAndPreserve(atoms);
}

namespace simulation
{
  /**
   * output information about the topology.
   */
  inline std::ostream & operator<<(std::ostream &os, topology &topo)
  {
    os << "a topology";
    return os;
  }
}

