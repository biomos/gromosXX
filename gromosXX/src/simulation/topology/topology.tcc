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

/**
 * @name the bond class
 * @{
 */
/**
 * get an iterator over the bonds.
 */
inline simulation::topology::bond::iterator simulation::topology::bond::begin()
{
  return iterator(m_bond_information);
}
/**
 * add a bond.
 */
inline void simulation::topology::bond::add(int i, int j, int type)
{
  bond_struct s;
  s.i = i;
  s.j = j;
  s.type = type;
  m_bond_information.push_back(s);
}

/**
 * iterator constructor.
 */
inline simulation::topology::bond::iterator
::iterator(std::vector<bond_struct> &bi)
{
  m_bond_it = bi.begin();
  m_bond_end = bi.end();
}

/**
 * end of list?
 */
bool inline simulation::topology::bond::iterator::eol()
{
  return m_bond_it == m_bond_end;
}

/**
 * increment.
 */
void inline simulation::topology::bond::iterator::operator++()
{
  ++m_bond_it;
}

/**
 * bond atom i.
 */
int inline simulation::topology::bond::iterator::i()
{
  return m_bond_it->i;
}

/**
 * bond atom j.
 */
int inline simulation::topology::bond::iterator::j()
{
  return m_bond_it->j;
}

/**
 * bond type.
 */
int inline simulation::topology::bond::iterator::type()
{
  return m_bond_it->type;
}

/**
 * @}
 */

/**
 * bond accessor
 */
inline simulation::topology::bond & simulation::topology::bonds()
{
  return m_bonds;
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

