/**
 * @file bond.tcc
 * inline methods for bond-topology.
 */

/**
 * get an iterator over the bonds.
 */
inline simulation::bond::iterator simulation::bond::begin()
{
  return iterator(m_bond_information);
}

/**
 * add a bond.
 */
inline void simulation::bond::add(int i, int j, int type)
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
inline simulation::bond::iterator
::iterator(std::vector<bond_struct> &bi)
{
  m_bond_it = bi.begin();
  m_bond_end = bi.end();
}

/**
 * end of list?
 */
bool inline simulation::bond::iterator::eol()
{
  return m_bond_it == m_bond_end;
}

/**
 * not end of list?
 */
bool inline simulation::bond::iterator::neol()
{
  return m_bond_it != m_bond_end;
}

/**
 * increment.
 */
void inline simulation::bond::iterator::operator++()
{
  ++m_bond_it;
}

/**
 * bond atom i.
 */
int inline simulation::bond::iterator::i()
{
  return m_bond_it->i;
}

/**
 * bond atom j.
 */
int inline simulation::bond::iterator::j()
{
  return m_bond_it->j;
}

/**
 * bond type.
 */
int inline simulation::bond::iterator::type()
{
  return m_bond_it->type;
}
