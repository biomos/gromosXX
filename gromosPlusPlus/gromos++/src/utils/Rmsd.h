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

// utils_Rmsd.h

#ifndef INCLUDED_UTILS_RMSD
#define INCLUDED_UTILS_RMSD


namespace gcore{
  class System;
}

namespace fit{
  class Reference;
  
}

namespace utils{
   class Property;
   class PropertyContainer;
  /**
   * Class Rmsd
   * A class to calculate the rmsd compared to a reference
   *
   * The atom positional root mean square deviation of the system 
   * is calculated with respect to a reference
   *
   * @class Rmsd
   * @author R. Buergi
   * @ingroup utils
   * @sa fit::Reference
   */
  class Rmsd{
    const fit::Reference *d_ref;
    const utils::PropertyContainer *d_prop_sys;
    const utils::PropertyContainer *d_prop_ref;
    // not implemented
    Rmsd();
    Rmsd(const Rmsd &);
    Rmsd &operator=(const Rmsd&);
  public:
    /**
     * Rmsd constructor
     */
    Rmsd(const fit::Reference *);    
    /**
     * Rmsd deconstructor
     */
    ~Rmsd(){}
    /** 
     * Method that actually calculate the positional rmsd between the reference
     * and the system
     * @param & System to calculate the rmsd for
     * @return the rmsd value
     */
    double rmsd(const gcore::System &);
    /** 
     * Method that actually calculate a property rmsd of a PropertyContainer
     * between the reference and the system
     * @param & System to calculate the rmsd for
     * @return the rmsd value
     */
    double rmsdproperty(const gcore::System &);

    /** 
     * Method to add a PropertyContainer for the rmsd calculation
     * @param * PropertyContainer of the reference system
     * @param * PropertyContainer of the system
     */
    void addproperty(const utils::PropertyContainer *, const utils::PropertyContainer *);

    typedef double (Rmsd::*MemPtr)(const gcore::System &);

  };
}

#endif
