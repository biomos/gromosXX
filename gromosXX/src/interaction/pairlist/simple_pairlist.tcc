/**
 * @file simple_pairlist.tcc
 * simple pairlist implementation.
 */

/**
 * clear all pairlists.
 */
template<typename t_simulation>
inline void interaction::simple_pairlist<t_simulation>
::clear_pairlist()
{
  m_pairlist.clear();
}

/**
 * Build up a pairlist, here
 * every particle interacts with all other particles.
 * To not count double, the pairlist contains only
 * interactions with particles with higher sequence
 * number.
 */
template<typename t_simulation>
void interaction::simple_pairlist<t_simulation>
::make_pairlist(t_simulation &simu)
{
  clear_pairlist();
  size_t num_atoms = simu.topology().num_solute_atoms();
  m_pairlist.resize(num_atoms);
  
  // now the a bit silly algorithm
  for(size_t i=0; i<num_atoms; ++i)
    for(size_t j=i+1; j<num_atoms; ++j)
      m_pairlist[i].push_back(j);
  
}

template<typename t_simulation>
void interaction::simple_pairlist<t_simulation>
::print_pairlist(std::ostream &os)
{
  os << "pairlist\n";
  os << "--------\n";
  
  std::vector< std::vector<int> >::const_iterator
    i = m_pairlist.begin(),
    i_to = m_pairlist.end();
  
  std::vector<int>::const_iterator j, j_to;
  
  for(int n=0; i != i_to; ++i, ++n){
    j = i->begin();
    j_to = i->end();
    os << std::setw(6) << n << " | ";
    for( ; j != j_to; ++j){
      os << std::setw(6) << *j;
    }
    os << std::endl;
  }
}

template<typename t_simulation>
inline typename interaction::simple_pairlist<t_simulation>::iterator 
interaction::simple_pairlist<t_simulation>::begin()
{
  return iterator(m_pairlist);
}

template<typename t_simulation>
inline interaction::simple_pairlist<t_simulation>::iterator::iterator
(std::vector< std::vector<int> > &pl)
  : m_pairlist(pl),
    m_atom(0)
{
  m_i = m_pairlist.begin();
  if (m_i != m_pairlist.end())
    m_j = m_i->begin();
}

template<typename t_simulation>
inline void interaction::simple_pairlist<t_simulation>::iterator::operator++()
{
  ++m_j;
  while (m_j == m_i->end()){
    ++m_i;
    ++m_atom;
    if (m_i == m_pairlist.end()) break;
    m_j = m_i->begin();
  }
}

template<typename t_simulation>
inline int interaction::simple_pairlist<t_simulation>::iterator::i()
{
  return m_atom;
}

template<typename t_simulation>
inline int interaction::simple_pairlist<t_simulation>::iterator::j()
{
  return *m_j;
}

template<typename t_simulation>
inline bool interaction::simple_pairlist<t_simulation>::iterator::eol()
{
  return m_i == m_pairlist.end();
}

  
  

								      



