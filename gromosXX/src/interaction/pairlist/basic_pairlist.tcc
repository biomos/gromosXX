/**
 * @file basic_pairlist.tcc
 * template methods of Basic_Pairlist
 */

template<typename t_simulation>
inline
interaction::Basic_Pairlist<t_simulation>
::iterator::iterator(basic_pairlist_type &pl)
  :   m_i(pl.begin()),
      m_j(m_i->begin()),
      m_pairlist(pl)
{
  
  // go to the first pair
  while(m_i != pl.end()){
    if (m_j == m_i->end()){
      ++m_i;
      m_j = m_i->begin();
    }
    else break;
  }
}

template<typename t_simulation>
inline void
interaction::Basic_Pairlist<t_simulation>::iterator::operator++()
{
  if(++m_j == m_i->end()){
    while(m_i !=  m_pairlist.end()){
      ++m_i;
      m_j = m_i->begin();
      if (m_j != m_i->end())
	break;
    }
  }
}

template<typename t_simulation>
inline bool
interaction::Basic_Pairlist<t_simulation>::iterator
::operator==(Basic_Pairlist<t_simulation>::iterator &it)
{
  if (m_i == it.m_i) return true;
  return false;
}

template<typename t_simulation>
inline bool
interaction::Basic_Pairlist<t_simulation>::iterator
::operator!=(typename Basic_Pairlist<t_simulation>::iterator &it)
{
  if (m_i != it.m_i) return true;
  return false;
}

template<typename t_simulation>
inline void
interaction::Basic_Pairlist<t_simulation>::iterator::row(unsigned int i)
{
  m_i = m_pairlist.begin() + i;
  if (m_i != m_pairlist.end())
    m_j = m_i->begin();
}

template<typename t_simulation>
inline typename interaction::Basic_Pairlist<t_simulation>::iterator
interaction::Basic_Pairlist<t_simulation>::begin()
{
  return iterator(*this);
}

template<typename t_simulation>
inline typename interaction::Basic_Pairlist<t_simulation>::iterator
interaction::Basic_Pairlist<t_simulation>::end()
{
  iterator it(*this);
  it.row(size());
  return it;
}

