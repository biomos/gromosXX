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

// utils_Neighbours.h

// Class that contains bonded neighbours of atom i

#ifndef INCLUDED_UTILS_NEIGHBOURS
#define INCLUDED_UTILS_NEIGHBOURS

#include <vector>

namespace gcore{
  class Simulation;
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace utils{
  /**
   * Class Neighbours
   * A vector of integers that contains all the bonded neighbours of an
   * atom
   *
   * Specifying the atom in different ways, this class is a vector of its
   * bonded immediate neighbours. Access is like for any vector of integers
   *
   * @class Neighbours
   * @author R. Buergi
   * @ingroup utils
   */
class Neighbours: public std::vector<int>{
  // not implemented
  Neighbours (const Neighbours&);
  Neighbours();
  Neighbours &operator=(const Neighbours&);
 public:
  /**
   * Neighbours constructor
   * @param sys The complete system is specified
   * @param mol The molecule number of the atom you are interested in
   * @param i   The atom number of the atom you are interested in
   */
  Neighbours(const gcore::System &sys, int mol, int i);
  /**
   * Neighbours constructor, expecially for solvent. 
   * @author M.A. Kastenholz
   * @param sys The complete system is specified
   * @param mol The number of the solvent. Should in this case always be 0
   * @param i   The number of the atom you are interested in
   * @param j   Any number this is just an ugly way to distinguish a solvent
   *            Neighbour construction from an solute. Don't ask why.
   */
  Neighbours(const gcore::System &sys, int mol, int i, int j);
  /**
   * Neighbours constructor
   * @param mol In this case you specify a single Molecule
   * @param k   And the number of the atom you are interested in
   */
  Neighbours(const gcore::Molecule &mol, int k);
  /**
   * Neighbours constructor
   * @param mol Specifying a MoleculeTopology is also possible
   * @param k   And the number of the atom you are interested in
   */
  Neighbours(const gcore::MoleculeTopology &mol, int k);
};
}
#endif
