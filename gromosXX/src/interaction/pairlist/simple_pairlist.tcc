/**
 * @file simple_pairlist.tcc
 * simple pairlist implementation.
 */

#ifndef INCLUDED_STDEXCEPT
#include <stdexcept>
#define INCLUDED_STDEXCEPT
#endif

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
void interaction::twin_range_pairlist<t_simulation>
::update(t_simulation &sim)
{
  short_range().clear();
  long_range().clear();

  size_t num_atoms = sim.topology().num_solute_atoms();

  short_range().resize(num_atoms);
  long_range().resize(num_atoms);

  math::VArray &pos = sim.system().pos();

  double d = 0;

  // square of the cutoff...
  double rcutp2 = sim.nonbonded_cutoff_short();
  rcutp2 *= rcutp2;
  
  double rcutl2 = sim.nonbonded_cutoff_long();
  rcutl2 *= rcutl2;

  math::Vec p;
  
  for(size_t i=0; i<num_atoms; ++i) {
    for(size_t j=i+1; j<num_atoms; ++j) {

      sim.system().periodicity().nearest_image(pos(i), pos(j), p);
      d = dot(p, p);

      if (d > rcutl2) 
        continue;
      else if (d > rcutp2)
        long_range()[i].push_back(j);
      else {
	// check it is not excluded
	if(i < sim.topology().solute().num_atoms())
	  if (sim.topology().all_exclusion(i).count(j))
	    continue;

        short_range()[i].push_back(j);
      }
    }
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
  os << "printing pairlist\n";
  os << "-----------------" << endl;

  typename simple_pairlist<t_simulation>::iterator it = pl.begin();
  while (it != pl.end()) {
    if (!it.j())
      os << endl << std::setw(6) << it.i() << " | " << flush;
    os << std::setw(6) << *it; 
    ++it;
  }
  os << std::endl;

  return os;
}
