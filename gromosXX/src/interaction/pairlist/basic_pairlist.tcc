/**
 * @file basic_pairlist.tcc
 * template methods of Basic_Pairlist
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE pairlist

#include "../../debug.h"

template<typename t_simulation, typename t_pairlist_algorithm>
inline
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::iterator(basic_pairlist_type &pl)
  :  
      m_pairlist(pl)
{

  DEBUG(7, "Pairlist size here " << pl.size());
  
  m_i = pl.begin();
  // go to the first pair
  while(m_i != pl.end()){
    m_j = m_i->begin();
    if (m_j == m_i->end()){
      ++m_i;
    }
    else {
      DEBUG(7, "we got a pair");
      break;
    }
  }
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline void
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::operator++()
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

template<typename t_simulation, typename t_pairlist_algorithm>
inline bool
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>::iterator
::operator==(typename Basic_Pairlist<t_simulation, t_pairlist_algorithm>
	     ::iterator &it)
{
  if (m_i == it.m_i) return true;
  return false;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline bool
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>::iterator
::operator!=(typename Basic_Pairlist<t_simulation, t_pairlist_algorithm>
	     ::iterator &it)
{
  if (m_i != it.m_i) return true;
  return false;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline void
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::row(unsigned int i)
{
  m_i = m_pairlist.begin() + i;
  DEBUG(7, "m_i set");
  if (m_i != m_pairlist.end()){
    DEBUG(7, "not the end");
    m_j = m_i->begin();
  }
  
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline unsigned int
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::i()
{
  return (m_i - m_pairlist.begin());
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline unsigned int
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::j()
{
  return *m_j;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline unsigned int
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::operator*()
{
  return *m_j;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::Basic_Pairlist(interaction::Nonbonded_Base &base)
  : t_pairlist_algorithm(*this, base)
{
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline typename interaction::Basic_Pairlist<t_simulation,
					    t_pairlist_algorithm>::iterator
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>::begin()
{
  DEBUG(7, "\tsize ot this : " << this->size());
  
  return iterator(*this);
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline typename interaction::Basic_Pairlist<t_simulation,
					    t_pairlist_algorithm>::iterator
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>::end()
{
  DEBUG(7, "pairlist::end");
  iterator it(*this);
  it.row(size());
  return it;
}

namespace interaction
{
  template<typename t_simulation, typename t_pairlist_algorithm>
  std::ostream & 
  operator<<(std::ostream &os, Basic_Pairlist<t_simulation, 
	     t_pairlist_algorithm> &pl)
  {
    os << "Printing pairlist" << std::endl;
    
    typename Basic_Pairlist<t_simulation, t_pairlist_algorithm>::iterator
      it = pl.begin(),
      to = pl.end();

    int ind = 1;
    while (it != to) {
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
  
}

