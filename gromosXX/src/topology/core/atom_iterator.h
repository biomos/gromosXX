/**
 * @file atom_iterator.h
 * iterator over atoms (of a subgroup)
 */

#ifndef INCLUDED_ATOM_ITERATOR_H
#define INCLUDED_ATOM_ITERATOR_H

namespace topology
{
  /**
   * @class Atom_Iterator
   * iterator over (a group of) atoms.
   */
  class Atom_Iterator
  {
  public:
    Atom_Iterator(size_t const atom)
      : m_atom(atom)
    {
    }
    bool operator==(Atom_Iterator const &it)const
    {
      return m_atom == it.m_atom;
    }
    bool operator!=(Atom_Iterator const &it)const
    {
      return m_atom != it.m_atom;
    }
    void operator++()
    {
      ++m_atom;
    }
    Atom_Iterator operator+(size_t n)const
    {
      Atom_Iterator dummy(*this);
      dummy.m_atom += n;
      return dummy;
    }
    
    size_t operator*()const
    {
      return m_atom;
    }
  private:
    size_t m_atom;
  };
    
} // topology

#endif