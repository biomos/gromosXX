/**
 * @file atom_group_iterator.h
 * iterator over a group of atoms.
 */

#ifndef INCLUDED_ATOM_GROUP_ITERATOR
#define INCLUDED_ATOM_GROUP_ITERATOR

namespace simulation
{
  /**
   * @class Atom_Group_Iterator
   * iterator over a group of atoms.
   */
  template<typename t_group_it>
  class Atom_Group_Iterator
  {
  public:
    /**
     * Constructor.
     */
    Atom_Group_Iterator(t_group_it const it)
      : m_it(it)
    {
    }
  
    bool operator==(Atom_Group_Iterator const &it)const
    {
      return m_it == it.m_it;
    }
    bool operator!=(Atom_Group_Iterator const &it)const
    {
      return !(*this == it);
    }
    void operator++()
    {
      ++m_it;
    }
    void operator+=(size_t const n)
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
    size_t num_atoms()const
    {
      return *(m_it+1) - *m_it;
    }
  
  private:
    t_group_it m_it;
  };
 
} // simulation

#endif
