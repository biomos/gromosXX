/**
 * @file simple_pairlist.tcc
 * simple pairlist implementation.
 */

template<typename t_simulation, typename t_filter>
inline 
interaction::simple_pairlist<t_simulation, t_filter>::iterator::iterator(
  t_pl_matrix &pl,
  int ii,
  int jj
)
  : 
  m_pairlist(pl),
  m_i(ii),
  m_j(jj)
{}

template<typename t_simulation, typename t_filter>
inline 
void 
interaction::simple_pairlist<t_simulation, t_filter>::iterator::operator++()
{
  ++j();
  if (j() == m_pairlist[i()].size()) {
    j() = 0;
    do {
      ++i();
      if (i() == m_pairlist.size()) 
        break;
    } while (!m_pairlist[i()].size());
  }
}

template<typename t_simulation, typename t_filter>
inline bool
interaction::simple_pairlist<t_simulation, t_filter>::iterator
::operator==(typename simple_pairlist<t_simulation, t_filter>::iterator it)
{
  return (
    &pairlist() == &(it.pairlist()) &&
    i() == it.i() &&
    j() == it.j()
  );
}

template<typename t_simulation, typename t_filter>
inline bool
interaction::simple_pairlist<t_simulation>::iterator
::operator!=(typename simple_pairlist<t_simulation, t_filter>::iterator it)
{
  return !operator==(it);
}

template<typename t_simulation, typename t_filter>
void interaction::simple_pairlist<t_simulation, t_filter>
::update(t_simulation &sim)
{

  clear();
  const size_t num_atoms = sim.topology().num_atoms();
  resize(num_atoms);

  const size_t num_solute_atoms = sim.topology().num_solute_atoms();
  
  // now the a bit silly algorithm
  size_t i, j;

  // solute
  for(i=0; i<num_solute_atoms; ++i){
    for(j=i+1; j<num_solute_atoms; ++j){
      // both are solute, check exclusions
      if (solute_pair(sim, i, j)) continue;
      
      (*this)[i].push_back(j);
    }
    for( ; j < num_atoms; ++j){
      // solute - solvent
      (*this)[i].push_back(j);
    }
  }
  // solvent
  for( ; i<num_atoms; ++i){
    for(j=i+1; j<num_atoms; ++j){
      // both are solvent, check whether same molecule
      if (solvent_pair(sim, i, j)) continue;
      
      (*this)[i].push_back(j);
    }
  }
    
}

template<typename t_simulation, typename t_filter>
inline 
typename interaction::simple_pairlist<t_simulation, t_filter>::iterator 
interaction::simple_pairlist<t_simulation, t_filter>::begin()
{
  for (unsigned int ii = 0; ii < size(); ii++)
    if ((*this)[ii].size())
      return iterator(*this, ii);
  return end();
}

template<typename t_simulation, typename t_filter>
inline 
typename interaction::simple_pairlist<t_simulation, t_filter>::iterator 
interaction::simple_pairlist<t_simulation, t_filter>::end()
{
  return iterator(*this, size());
}

template<typename t_simulation, typename t_filter>
std::ostream& 
interaction::operator<<(
  std::ostream &os, 
  class simple_pairlist<t_simulation, t_filter>& pl
)
{
  // os << "printing pairlist\n";
  // os << "-----------------" << endl;

  typename simple_pairlist<t_simulation, t_filter>::iterator it = pl.begin();
  int ind = 1;
  while (it != pl.end()) {
    if (!it.j()){
      ind = 1;
      os << endl << std::setw(5) << it.i() << ": " << flush;
    }
    os << std::setw(5) << *it << " "; 
    if (!(ind%15)) os << "\n\t";
    ++it;
    ++ind;
  }
  os << std::endl;

  return os;
}
