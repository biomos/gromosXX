/**
 * @file simple_pairlist.tcc
 * simple pairlist implementation.
 */

#ifndef INCLUDED_STDEXCEPT
#include <stdexcept>
#define INCLUDED_STDEXCEPT
#endif

#ifndef INCLUDED_MATH_PERIODICITY
#include "math/periodicity.h"
#define INCLUDED_MATH_PERIODICITY
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
::update(t_simulation &simu)
{

  clear();
  size_t num_atoms = simu.topology().num_solute_atoms();
  resize(num_atoms);
  
  // now the a bit silly algorithm
  for(size_t i=0; i<num_atoms; ++i)
    for(size_t j=i+1; j<num_atoms; ++j)
      (*this)[i].push_back(j);
}

template<typename t_simulation>
void interaction::twin_range_pairlist<t_simulation>
::update(t_simulation &simu)
{
  short_range().clear();
  long_range().clear();

  size_t num_atoms = simu.topology().num_solute_atoms();

  short_range().resize(num_atoms);
  long_range().resize(num_atoms);

  math::VArray p = simu.system().pos();

  double d = 0;

  double rcutp = 4;
  double rcutl = 5;

  for(size_t i=0; i<num_atoms; ++i) {
    for(size_t j=i+1; j<num_atoms; ++j) {
      d = sqrt(sum(p(j) - p(i)));
      if (d > rcutl) 
        continue;
      else if (d > rcutp)
        long_range()[i].push_back(j);
      else 
        short_range()[i].push_back(j);
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
  os << "--------" << endl;

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
