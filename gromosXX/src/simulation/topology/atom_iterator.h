/**
 * @file atom_iterator.h
 * iterator over atoms (of a subgroup)
 */

#ifndef INCLUDED_ATOM_ITERATOR_H
#define INCLUDED_ATOM_ITERATOR_H

namespace simulation
{
  /**
   * @class Atom_Iterator
   * iterator over (a group of) atoms.
   */
  class Atom_Iterator
  {
  public:
    Atom_Iterator(size_t atom)
      : m_atom(atom)
    {
    }
    bool operator==(Atom_Iterator &it)
    {
      return m_atom == it.m_atom;
    }
    bool operator!=(Atom_Iterator &it)
    {
      return m_atom != it.m_atom;
    }
    void operator++()
    {
      ++m_atom;
    }
    size_t operator*()
    {
      return m_atom;
    }
  private:
    size_t m_atom;
  };
    
} // simulation

#endif
