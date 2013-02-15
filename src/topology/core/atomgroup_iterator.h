/**
 * @file atomgroup_iterator.h
 * iterator over a group of atoms.
 */

#ifndef INCLUDED_ATOMGROUP_ITERATOR
#define INCLUDED_ATOMGROUP_ITERATOR

namespace topology
{
  /**
   * @class Atomgroup_Iterator
   * iterator over a group of atoms.
   */
  template<typename t_group_it>
  class Atomgroup_Iterator
  {
  public:
    /**
     * Constructor.
     */
    Atomgroup_Iterator(t_group_it const it)
      : m_it(it)
    {
    }

    Atomgroup_Iterator()
    {
    }
    
    bool operator==(Atomgroup_Iterator const &it)const
    {
      return m_it == it.m_it;
    }
    bool operator<(Atomgroup_Iterator const &it)const
    {
      return m_it < it.m_it;
    }
    bool operator!=(Atomgroup_Iterator const &it)const
    {
      return !(*this == it);
    }
    void operator++()
    {
      ++m_it;
    }
    void operator+=(unsigned int n)
    {
      m_it += n;
    }
    
    t_group_it const & operator*()const
    {
      return m_it;
    }
    /**
     * maybe move this one out?
     * is it already too specified?
     */
    Atom_Iterator begin()const
    {
      return Atom_Iterator(*m_it);
    }
    Atom_Iterator end()const
    {
      return Atom_Iterator(*(m_it+1));
    }
    unsigned int num_atoms()const
    {
      return *(m_it+1) - *m_it;
    }
  
  private:
    t_group_it m_it;
  };
 
} // topology

#endif
