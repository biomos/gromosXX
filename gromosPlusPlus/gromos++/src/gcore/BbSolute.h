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
#ifndef INCLUDED_BBSOLUTE
#define INCLUDED_BBSOLUTE

#include <string>
#include <vector>

#include "Exclusion.h"
#include "MoleculeTopology.h"

namespace gcore{

  class GromosForceField;
  class AtomTopology;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  /**
   * Class BbSolute
   * Purpose: defines a Building block for solute molecules (MTBUILDBLSOLUTE)
   *
   * Description:
   * Defines a Building block that is (part of) a solute molecule. Note 
   * that for the atoms there is a direct accessor, but that the bonds,
   * angles etc. are accessed via the iterators, just as in the 
   * MoleculeTopology.
   *
   * @class BbSolute
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::MoleculeTopology
   * @sa gcore::BuildingBlock
   */
  class BbSolute: public gcore::MoleculeTopology{
    std::vector<Exclusion> d_pexcl;
    int d_rep;
    
  public:
    /**
     * BbSolute Constructor
     */
    BbSolute(){d_rep=0;};
    /**
     * BbSolute copy constructor
     */
    BbSolute(const BbSolute &);
    /**
     * MoleculeTopology copy constructor
     */
    BbSolute(const gcore::MoleculeTopology &);
    
    /**
     * BbSolute deconstructor
     */
    ~BbSolute(){};
    /**
     * Member operator = to copy one BbSolute into the other
     */
    BbSolute &operator=(const BbSolute &);

    // Methods
    /**
     * Member function to add previous exclusions to the building block.
     * Previous exclusions are the exclusions that should be applied to the
     * previous building block
     * @param a Exclusions of type gcore::Exclusion
     */
    void addPexcl(const Exclusion &a);
    /**
     * Member function to set the name of the building block
     * @param s String containing the name
     */
    void setResName(const std::string &s);
    /**
     * Member function to set the number of atoms that replace atoms in
     * linked building blocks. If i &lt; 0 the first i atoms replace the last
     * i atoms of the previous building block. If i &gt; 0 the last i atoms
     * replace the first i atoms of the following building block.<br>
     * @see make_top.
     */
    void setRep(const int i);

    // Accessors
    /** 
     * Accessor, returns the number of preceding exclusions in the building
     * block
     */
    int numPexcl()const;
    /** 
     * Accessor, returns the i-th set of preceding exclusions
     */
    const Exclusion& pexcl(int i) const;
    /**
     * Accessor, returns the name of the building block
     */
    const std::string &resName()const;
    /**
     * Accessor, returns the number of atoms to be replaced. If i &lt; 0 the
     * first i atoms replace the last i atoms of the previous building
     * block. If i &gt; 0 the last i atoms replace the first i atoms of the
     * following building block.
     * @see make_top.
     */
    int rep()const;

  }; /* class BbSolute */
} /* Namespace */ 



#endif
