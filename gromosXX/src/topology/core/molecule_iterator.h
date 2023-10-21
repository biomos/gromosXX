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
 * @file molecule_iterator.h
 * iterator over the molecules.
 */

#ifndef INCLUDED_MOLECULE_ITERATOR
#define INCLUDED_MOLECULE_ITERATOR

namespace topology
{
  /**
   * @class Molecule_Iterator
   * iterator over the molecules (SUBMOLECULES block)
   */
  class Molecule_Iterator 
    : public Atomgroup_Iterator<std::vector<unsigned int>::const_iterator>
  {
  public:
    Molecule_Iterator(std::vector<unsigned int>::const_iterator mol_it)
      : Atomgroup_Iterator<std::vector<unsigned int>::const_iterator>(mol_it) {};
    

  private:
  };
  
} // topology

#endif

  
