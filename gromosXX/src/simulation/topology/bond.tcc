/**
 * @file bond.tcc
 * inline methods for bond-topology.
 */

/**
 * get an iterator over the bonds.
 */
inline simulation::Bond::iterator simulation::Bond::begin()
{
  return iterator(m_bond_information);
}

/**
 * add a bond.
 */
inline void simulation::Bond::add(int i, int j, int type)
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
inline simulation::Bond::iterator
::iterator(std::vector<bond_struct> &bi)
{
  m_bond_it = bi.begin();
  m_bond_end = bi.end();
}

/**
 * end of list?
 */
bool inline simulation::Bond::iterator::eol()
{
  return m_bond_it == m_bond_end;
}

/**
 * not end of list?
 */
bool inline simulation::Bond::iterator::neol()
{
  return m_bond_it != m_bond_end;
}

/**
 * increment.
 */
void inline simulation::Bond::iterator::operator++()
{
  ++m_bond_it;
}

/**
 * bond atom i.
 */
int inline simulation::Bond::iterator::i()
{
  return m_bond_it->i;
}

/**
 * bond atom j.
 */
int inline simulation::Bond::iterator::j()
{
  return m_bond_it->j;
}

/**
 * bond type.
 */
int inline simulation::Bond::iterator::type()
{
  return m_bond_it->type;
}
