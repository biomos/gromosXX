/**
 * @file pairlist.tcc
 * inline methods of Pairlist
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE pairlist

#include "../../debug.h"

inline
interaction::Pairlist
::iterator::iterator(std::vector<std::vector<size_t> > &pl)
  : m_pairlist(pl)
{

  DEBUG(7, "Pairlist size @iteractor construction " << pl.size());
  
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

inline void
interaction::Pairlist::iterator
::operator++()
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

inline bool
interaction::Pairlist::iterator
::operator==(Pairlist::iterator &it)
{
  if (m_i == it.m_i) return true;
  return false;
}

inline bool
interaction::Pairlist::iterator
::operator!=(Pairlist::iterator &it)
{
  if (m_i != it.m_i){
    return true;
  }
  return false;
}

inline void
interaction::Pairlist::iterator
::row(unsigned int i)
{
  DEBUG(7, "pairlist iterator row() " << i);

  m_i = m_pairlist.begin() + i;
  if (m_i != m_pairlist.end()){
    m_j = m_i->begin();
  }
}

inline unsigned int
interaction::Pairlist::iterator
::i()
{
  DEBUG(10, "Pairlist::i()");
  return (m_i - m_pairlist.begin());
}

inline unsigned int
interaction::Pairlist::iterator
::j()
{
  DEBUG(10, "Pairlist::j()");
  return *m_j;
}

inline unsigned int
interaction::Pairlist::iterator
::operator*()
{
  DEBUG(10, "Pairlist::operator*()");
  return *m_j;
}

inline
interaction::Pairlist
::Pairlist()
{
}

inline interaction::Pairlist::iterator
interaction::Pairlist::begin()
{
  return iterator(*this);
}

inline interaction::Pairlist::iterator
interaction::Pairlist::end()
{
  iterator it(*this);
  it.row(size());
  return it;
}

namespace interaction
{
  std::ostream & 
  operator<<(std::ostream &os, Pairlist &pl)
  {
    {
      os << "Pairlist" << std::endl;
    
      Pairlist::iterator
	it = pl.begin(),
	to = pl.end();

      size_t ii = 999999999;

      int ind = 1;
      while (it != to) {
	if (ii != it.i()){
	  ii = it.i();
	  ind = 1;
	  os << std::endl << std::setw(5) << it.i() << ": " << std::flush;
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
  
} // namespace interaction

