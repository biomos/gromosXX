/**
 * @file improper_dihedral.tcc
 * inline methods for improper_dihedral-topology.
 */

/**
 * get an iterator over the improper_dihedrals.
 */
inline simulation::Improper_dihedral::iterator 
simulation::Improper_dihedral::begin()
{
  return iterator(m_improper_dihedral_information);
}

/**
 * add an improper dihedral.
 */
inline void simulation::Improper_dihedral::add(int i, int j, int k, int l, 
					       int type)
{
  improper_dihedral_struct s;
  s.i = i;
  s.j = j;
  s.k = k;
  s.l = l;
  s.type = type;
  m_improper_dihedral_information.push_back(s);
}

/**
 * iterator constructor.
 */
inline simulation::Improper_dihedral::iterator
::iterator(std::vector<improper_dihedral_struct> &ii)
{
  m_improper_dihedral_it = ii.begin();
  m_improper_dihedral_end = ii.end();
}

/**
 * end of list?
 */
bool inline simulation::Improper_dihedral::iterator::eol()
{
  return m_improper_dihedral_it == m_improper_dihedral_end;
}

/**
 * not end of list?
 */
bool inline simulation::Improper_dihedral::iterator::neol()
{
  return m_improper_dihedral_it != m_improper_dihedral_end;
}

/**
 * increment.
 */
void inline simulation::Improper_dihedral::iterator::operator++()
{
  ++m_improper_dihedral_it;
}

/**
 * improper dihedral atom i.
 */
int inline simulation::Improper_dihedral::iterator::i()
{
  return m_improper_dihedral_it->i;
}

/**
 * improper dihedral atom j.
 */
int inline simulation::Improper_dihedral::iterator::j()
{
  return m_improper_dihedral_it->j;
}
/**
 * improper dihedral atom k.
 */
int inline simulation::Improper_dihedral::iterator::k()
{
  return m_improper_dihedral_it->k;
}
/**
 * improper dihedral atom l.
 */
int inline simulation::Improper_dihedral::iterator::l()
{
  return m_improper_dihedral_it->l;
}

/**
 * improper dihedral type.
 */
int inline simulation::Improper_dihedral::iterator::type()
{
  return m_improper_dihedral_it->type;
}
