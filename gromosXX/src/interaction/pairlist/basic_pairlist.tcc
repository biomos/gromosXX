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
  DEBUG(7, "pairlist iterator constructor: end == " << int(m_i == pl.end()));
  
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
  if (m_i != it.m_i){
    DEBUG(10, "pairlist iterators are !=");
    return true;
  }
  DEBUG(10, "pairlist iterators are ==");  
  return false;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline void
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::row(unsigned int i)
{
  DEBUG(7, "pairlist iterator row() " << i);

  m_i = m_pairlist.begin() + i;
  DEBUG(7, "m_i set");
  if (m_i != m_pairlist.end()){
    DEBUG(7, "not the end");
    m_j = m_i->begin();
  }
  DEBUG(7, "the end");
  
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline unsigned int
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::i()
{
  DEBUG(10, "Basic_Pairlist::i()");
  return (m_i - m_pairlist.begin());
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline unsigned int
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::j()
{
  DEBUG(10, "Basic_Pairlist::j()");
  return *m_j;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline unsigned int
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::iterator::operator*()
{
  DEBUG(10, "Basic_Pairlist::operator*()");
  return *m_j;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::Basic_Pairlist(interaction::Nonbonded_Base &base)
  : t_pairlist_algorithm(*this, m_perturbed_pairlist, base)
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

template<typename t_simulation, typename t_pairlist_algorithm>
inline typename 
interaction::Basic_Pairlist<t_simulation,
			    t_pairlist_algorithm>::iterator
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::perturbed_begin()
{
  DEBUG(7, "\tsize ot perturbed pairlist : " 
	<< m_perturbed_pairlist.size());
  
  return iterator(m_perturbed_pairlist);
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline typename interaction::Basic_Pairlist<t_simulation,
					    t_pairlist_algorithm>::iterator
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::perturbed_end()
{
  iterator it(m_perturbed_pairlist);
  it.row(m_perturbed_pairlist.size());
  return it;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline interaction::basic_pairlist_type &
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::perturbed()
{
  return m_perturbed_pairlist;
}

template<typename t_simulation, typename t_pairlist_algorithm>
inline interaction::basic_pairlist_type const &
interaction::Basic_Pairlist<t_simulation, t_pairlist_algorithm>
::perturbed()const
{
  return m_perturbed_pairlist;
}

namespace interaction
{
  template<typename t_simulation, typename t_pairlist_algorithm>
  std::ostream & 
  operator<<(std::ostream &os, Basic_Pairlist<t_simulation, 
	     t_pairlist_algorithm> &pl)
  {
    {
      
      os << "Basic_Pairlist" << std::endl;
    
      typename Basic_Pairlist<t_simulation, t_pairlist_algorithm>::iterator
	it = pl.begin(),
	to = pl.end();

      size_t ii = -1;

      int ind = 1;
      while (it != to) {
	if (ii != it.i()){
	  ii = it.i();
	  ind = 1;
	  os << endl << std::setw(5) << it.i() << ": " << flush;
	}
	os << std::setw(5) << *it << " "; 
	if (!(ind%15)) os << "\n\t";
	++it;
	++ind;
      }
      os << std::endl;
    }
    
    if (pl.perturbed().size()){
      os << "Basic_Pairlist: Perturbed atoms" << std::endl;
      
      typename Basic_Pairlist<t_simulation, t_pairlist_algorithm>::iterator
	it = pl.perturbed_begin(),
	to = pl.perturbed_end();

      size_t ii = -1;
      int ind = 1;
      while (it != to) {
	if (ii != it.i()){
	  ii = it.i();
	  ind = 1;
	  os << endl << std::setw(5) << it.i() << ": " << flush;
	}
	os << std::setw(5) << *it << " "; 
	if (!(ind%15)) os << "\n\t";
	++it;
	++ind;
      }
      os << std::endl;
    }
    
    return os;
  }
  
}

