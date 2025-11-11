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



// utils_FfExpert.h

// Class that contains statistically information on the occurrence
// of types and charges in a force field building block

#ifndef INCLUDED_UTILS_FFEXPERT
#define INCLUDED_UTILS_FFEXPERT


#include <map>
#include <string>
#include <vector>

#include "FfExpertGraph.h"
#include "../gmath/Vec.h"
#include "../gcore/Angle.h"
#include "../gcore/BuildingBlock.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/BbSolute.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Dihedral.h"
#include "../gcore/Improper.h"
#include "../gcore/MoleculeTopology.h"

namespace gcore
{
  class BuildingBlock;
  class Bond;
  class Angle;
  class Improper;
  class Dihedral;
}

namespace utils
{
  /**
   * Class FfExpert
   * contains statistical data on the occurence of types in a BuildingBlock
   * usefull for checking of consistency and suggesting types.
   *
   * Description:
   * This is a low-level expert system that knows which types of bonds, 
   * angles etc. are most commonly seen with certain atomic Iac values. The 
   * IAC are connected to the first letters of atom names
   *
   * @class FfExpert
   * @author C. Oostenbrink
   * @ingroup utils
   * @sa gcore::BuildingBlock
   */
  class FfExpert{
  public:
    struct counter
    {
      int type, occurence;
      counter(int i, int j)
      {
	type=i;
	occurence = j;
      }
      
    };
  protected:
    std::multimap<std::string,counter> d_name2iac;
    std::multimap<int,counter> d_iac2mass;
    std::multimap<int,counter> d_iac2charge;
    std::multimap<gcore::Bond,counter> d_iac2bond;
    std::multimap<gcore::Angle,counter> d_iac2angle;
    std::multimap<gcore::Improper,counter> d_iac2improper;
    std::multimap<gcore::Dihedral,counter> d_iac2dihedral;
    std::vector<double> d_chargeType;
    std::vector<FfExpertGraph> d_graphs;
  public:
    /**
     * Standard constructor
     */
    FfExpert() = default;
    /**
     * Constructor with learning of a file
     */
    FfExpert(gcore::BuildingBlock const & mtb, const FfExpertGraphMapper * mapper = NULL) {
      learn(mtb, mapper);
    }
    /**
     * Function to learn about a BuildingBlock
     */
    void learn(gcore::BuildingBlock const & mtb, const FfExpertGraphMapper * mapper = NULL);
    /**
     * Accessor to the names
     */
    void name2iac(std::string s, std::vector<counter> &v) const ;
    /**
     * Accessor to the substructures
     */
    void substructure2iac(unsigned int i, const FfExpertGraph & query,
            std::vector<std::vector<utils::Vertex> > & iacs) const;
    /**
     * Accessor to the masses
     */
    void iac2mass(int i, std::vector<counter> &v) const;
    /**
     * Accessor to charge types
     */
    void iac2charge(int i, std::vector<counter> &v) const;
    /**
     * Accessor the the charge via the substructure
     */
    void substructure2charge(unsigned int i, const FfExpertGraph & query,
            std::vector<std::vector<utils::Vertex> > & charge);
    /**
     * Accessor to charge as function of type
     */
    double charge(int i) const;
    /**
     * Accessor to bonds
     */
    void iac2bond(gcore::Bond const & b, std::vector<counter> &v) const;
 /**
     * Accessor to angles
     */
    void iac2angle(gcore::Angle const & b, std::vector<counter> &v) const;
    /**
     * Accessor to impropers
     */
    void iac2improper(gcore::Improper const & b, std::vector<counter> &v) const;
 /**
     * Accessor to dihedrals
     */
    void iac2dihedral(gcore::Dihedral const & b, std::vector<counter> &v) const;
 

  };

int sort(std::vector<FfExpert::counter> &v, bool tt=true);

inline void FfExpert::name2iac(std::string s, std::vector<counter> &v) const {
  v.clear();
  if (d_name2iac.count(s)) {
    for (std::multimap<std::string, FfExpert::counter>::const_iterator
             iter = d_name2iac.lower_bound(s),
             to = d_name2iac.upper_bound(s);
         iter != to; ++iter) {
      v.push_back(iter->second);
    }
  }
}

inline void FfExpert::iac2mass(int i, std::vector<counter> &v) const {
  v.clear();
  if (d_iac2mass.count(i)) {
    for (std::multimap<int, FfExpert::counter>::const_iterator
             iter = d_iac2mass.lower_bound(i),
             to = d_iac2mass.upper_bound(i);
         iter != to; ++iter) {
      v.push_back(iter->second);
    }
  }
}

inline void FfExpert::iac2charge(int i, std::vector<counter> &v) const {
  v.clear();
  if (d_iac2charge.count(i)) {
    for (std::multimap<int, FfExpert::counter>::const_iterator
             iter = d_iac2charge.lower_bound(i),
             to = d_iac2charge.upper_bound(i);
         iter != to; ++iter) {
      v.push_back(iter->second);
    }
  }
}

inline double FfExpert::charge(int i) const {
  return d_chargeType[i];
}

inline void FfExpert::iac2bond(gcore::Bond const & b, std::vector<counter> &v) const {
  v.clear();
  if (d_iac2bond.count(b)) {
    for (std::multimap<gcore::Bond, FfExpert::counter>::const_iterator
             iter = d_iac2bond.lower_bound(b),
             to = d_iac2bond.upper_bound(b);
         iter != to; ++iter) {
      v.push_back(iter->second);
    }
  }
}

inline void FfExpert::iac2angle(gcore::Angle const &b,
                                std::vector<counter> &v) const {
  v.clear();
  if (d_iac2angle.count(b)) {
    for (std::multimap<gcore::Angle, FfExpert::counter>::const_iterator
             iter = d_iac2angle.lower_bound(b),
             to = d_iac2angle.upper_bound(b);
         iter != to; ++iter) {
      v.push_back(iter->second);
    }
  }
}

inline void FfExpert::iac2improper(gcore::Improper const &b,
                                   std::vector<counter> &v) const {
  v.clear();
  if (d_iac2improper.count(b)) {
    for (std::multimap<gcore::Improper, FfExpert::counter>::const_iterator
             iter = d_iac2improper.lower_bound(b),
             to = d_iac2improper.upper_bound(b);
         iter != to; ++iter) {
      v.push_back(iter->second);
    }
  }
}

inline void FfExpert::iac2dihedral(gcore::Dihedral const &b,
                                   std::vector<counter> &v) const {
  v.clear();
  if (d_iac2dihedral.count(b)) {
    for (std::multimap<gcore::Dihedral, FfExpert::counter>::const_iterator
             iter = d_iac2dihedral.lower_bound(b),
             to = d_iac2dihedral.upper_bound(b);
         iter != to; ++iter) {
      v.push_back(iter->second);
    }
  }
}
} // namespace utils

#endif
