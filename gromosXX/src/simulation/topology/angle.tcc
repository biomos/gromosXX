/**
 * @file angle.tcc
 * inline methods for bond-topology.
 */

/**
 * get an iterator over the bonds.
 */
inline simulation::Angle::iterator simulation::Angle::begin()
{
  return iterator(m_angle_information);
}

/**
 * add an angle.
 */
inline void simulation::Angle::add(int i, int j, int k, int type)
{
  angle_struct s;
  s.i = i;
  s.j = j;
  s.k = k;
  s.type = type;
  m_angle_information.push_back(s);
}

/**
 * iterator constructor.
 */
inline simulation::Angle::iterator
::iterator(std::vector<angle_struct> &bi)
{
  m_angle_it = bi.begin();
  m_angle_end = bi.end();
}

/**
 * end of list?
 */
bool inline simulation::Angle::iterator::eol()
{
  return m_angle_it == m_angle_end;
}

/**
 * not end of list?
 */
bool inline simulation::Angle::iterator::neol()
{
  return m_angle_it != m_angle_end;
}

/**
 * increment.
 */
void inline simulation::Angle::iterator::operator++()
{
  ++m_angle_it;
}

/**
 * angle atom i.
 */
int inline simulation::Angle::iterator::i()
{
  return m_angle_it->i;
}

/**
 * angle atom j.
 */
int inline simulation::Angle::iterator::j()
{
  return m_angle_it->j;
}
/**
 * angle atom k.
 */
int inline simulation::Angle::iterator::k()
{
  return m_angle_it->k;
}

/**
 * bond type.
 */
int inline simulation::Angle::iterator::type()
{
  return m_angle_it->type;
}
