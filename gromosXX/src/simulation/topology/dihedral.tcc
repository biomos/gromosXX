/**
 * @file dihedral.tcc
 * inline methods for dihedral-topology.
 */

/**
 * get an iterator over the dihedrals.
 */
inline simulation::Dihedral::iterator 
simulation::Dihedral::begin()
{
  return iterator(m_dihedral_information);
}

/**
 * add a dihedral.
 */
inline void simulation::Dihedral::add(int i, int j, int k, int l, 
				      int type)
{
  dihedral_struct s;
  s.i = i;
  s.j = j;
  s.k = k;
  s.l = l;
  s.type = type;
  m_dihedral_information.push_back(s);
}

/**
 * iterator constructor.
 */
inline simulation::Dihedral::iterator
::iterator(std::vector<dihedral_struct> &di)
{
  m_dihedral_it = di.begin();
  m_dihedral_end = di.end();
}

/**
 * end of list?
 */
bool inline simulation::Dihedral::iterator::eol()
{
  return m_dihedral_it == m_dihedral_end;
}

/**
 * not end of list?
 */
bool inline simulation::Dihedral::iterator::neol()
{
  return m_dihedral_it != m_dihedral_end;
}

/**
 * increment.
 */
void inline simulation::Dihedral::iterator::operator++()
{
  ++m_dihedral_it;
}

/**
 * dihedral atom i.
 */
int inline simulation::Dihedral::iterator::i()
{
  return m_dihedral_it->i;
}

/**
 * dihedral atom j.
 */
int inline simulation::Dihedral::iterator::j()
{
  return m_dihedral_it->j;
}
/**
 * dihedral atom k.
 */
int inline simulation::Dihedral::iterator::k()
{
  return m_dihedral_it->k;
}
/**
 * dihedral atom l.
 */
int inline simulation::Dihedral::iterator::l()
{
  return m_dihedral_it->l;
}

/**
 * improper dihedral type.
 */
int inline simulation::Dihedral::iterator::type()
{
  return m_dihedral_it->type;
}
