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

// gcore_BbSolute.h
#ifndef INCLUDED_BBLINK
#define INCLUDED_BBLINK

#include <vector>

#include "BbSolute.h"

namespace gcore{

  class GromosForceField;
  class AtomTopology;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  class BbSolute;
  /**
   * Class BbLINK
   * Purpose: defines a Building block for linking 
   *
   *
   * @class BbLink
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::MoleculeTopology
   * @sa gcore::BuildingBlock
   */
  class BbLink: public gcore::BbSolute{
    std::vector<int> d_linkres;
  public:
    /**
     * BbSolute Constructor
     */
    BbLink(){setRep(0);};
    /**
     * BbSolute copy constructor
     */
    BbLink(const BbLink &);
    /**
     * BbSolute deconstructor
     */
    ~BbLink(){};
    /**
     * Member operator = to copy one BbSolute into the other
     */
    BbLink &operator=(const BbLink &);

    // Methods
    /**
     * Member function to store the residue identifier for a given atom
     * @param a atom
     * @param i residue
     */
    void setLinkRes(const unsigned int a, const unsigned int i);

    // Accessors
    /** 
     * Accessor, returns the residue identifier for a given atom a
     */
    int linkRes(int a)const;

  }; /* class BbLink */
} /* Namespace */ 



#endif
