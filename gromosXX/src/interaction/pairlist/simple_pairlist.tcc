/**
 * @file simple_pairlist.tcc
 * simple pairlist implementation.
 */

template<typename t_simulation>
inline 
interaction::simple_pairlist<t_simulation>::iterator::iterator(
  t_pl_matrix &pl,
  int ii,
  int jj
)
  : 
  m_pairlist(pl),
  m_i(ii),
  m_j(jj)
{}

template<typename t_simulation>
inline 
void 
interaction::simple_pairlist<t_simulation>::iterator::operator++()
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

template<typename t_simulation>
inline 
bool
interaction::simple_pairlist<t_simulation>::iterator::operator==(
  typename simple_pairlist<t_simulation>::iterator it
)
{
  return (
    &pairlist() == &(it.pairlist()) &&
    i() == it.i() &&
    j() == it.j()
  );
}

template<typename t_simulation>
inline 
bool
interaction::simple_pairlist<t_simulation>::iterator::operator!=(
  typename simple_pairlist<t_simulation>::iterator it
)
{
  return !operator==(it);
}

template<typename t_simulation>
void interaction::simple_pairlist<t_simulation>
::update(t_simulation &sim)
{

  clear();
  size_t num_atoms = sim.topology().num_solute_atoms();
  resize(num_atoms);
  
  // now the a bit silly algorithm
  for(size_t i=0; i<num_atoms; ++i)
    for(size_t j=i+1; j<num_atoms; ++j){
      // check if not excluded
      if(i < sim.topology().solute().num_atoms())
	if (sim.topology().all_exclusion(i).count(j))
	  continue;

      (*this)[i].push_back(j);
    }
  
}

template<typename t_simulation>
inline 
typename interaction::simple_pairlist<t_simulation>::iterator 
interaction::simple_pairlist<t_simulation>::begin()
{
  for (unsigned int ii = 0; ii < size(); ii++)
    if ((*this)[ii].size())
      return iterator(*this, ii);
  return end();
}

template<typename t_simulation>
inline 
typename interaction::simple_pairlist<t_simulation>::iterator 
interaction::simple_pairlist<t_simulation>::end()
{
  return iterator(*this, size());
}

template<typename t_simulation>
std::ostream& 
interaction::operator<<(
  std::ostream &os, 
  class simple_pairlist<t_simulation>& pl
)
{
  // os << "printing pairlist\n";
  // os << "-----------------" << endl;

  typename simple_pairlist<t_simulation>::iterator it = pl.begin();
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
