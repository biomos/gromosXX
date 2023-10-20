/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
    Atom_Iterator(unsigned int atom)
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
    Atom_Iterator operator+(unsigned int n)const
    {
      Atom_Iterator dummy(*this);
      dummy.m_atom += n;
      return dummy;
    }
    
    unsigned int operator*()const
    {
      return m_atom;
    }
  private:
    unsigned int m_atom;
  };
    
} // topology

#endif
