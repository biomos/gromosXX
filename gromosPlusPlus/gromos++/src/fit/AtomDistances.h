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

// fit_AtomDistances.h

#ifndef INCLUDED_FIT_ATOMDISTANCES
#define INCLUDED_FIT_ATOMDISTANCES

#include <vector>

#include "../gromos/Exception.h"
#include "../gmath/Vec.h"

namespace fit{
  /**
   * Class AtomDistances
   * This class calculates the structural difference between to structures
   * based on the interatomic distances
   *
   * @class AtomDistances
   * @author Chris Oostenbrink
   * @ingroup fit
   * @sa fit::FastRotationalFit
   * @sa fit::TranslationalFit
   * @sa utils::Rmsd
   */
  class AtomDistances{
    /**
     * atoms to calculate the distances for
     */
    std::vector<bool> d_dist_spec;
    /**
     * number of atoms used in dist calculation
     */
    int d_dist_num_atoms;

  public:
    /**
     * Constructor: no specifications
     */
    AtomDistances() 
      : d_dist_spec(0), d_dist_num_atoms(0) {};
    /**
     * Constructor
     */
    AtomDistances(std::vector<bool> dist_spec)
      : d_dist_spec(dist_spec),
	d_dist_num_atoms(0){
      for(size_t i=0; i<d_dist_spec.size(); ++i)
	if(d_dist_spec[i]) ++d_dist_num_atoms;
    };
    
    /**
     * calculate the distance (using atoms specified for dist) on (reduced) atom positions
     * using the specified rotation matrix on sys.
     */
    double dist(std::vector<gmath::Vec> const &ref,
		std::vector<gmath::Vec> const &sys)const;
    

    /**
     * FastRotationalFit exception
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string &what): 
	gromos::Exception("RotationalFit",what){}
    };
    
  };

}

#endif    
