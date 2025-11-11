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

// fit_RotationalFit.h

#ifndef INCLUDED_FIT_ROTATIONALFIT
#define INCLUDED_FIT_ROTATIONALFIT


#include "../gromos/Exception.h"
#include "../utils/AtomSpecifier.h"

namespace gcore{
  class System;
}

namespace fit{
  class RotationalFit_i;
  class Reference;
  /**
   * Class RotationalFit
   * This class performs a rotational fit on the system
   *
   * A least squares fitting of one system is performed relative to a
   * reference system. The atoms that are taken into account are defined
   * by the Reference class
   *
   * @class RotationalFit
   * @author R. Buergi
   * @ingroup fit
   * @sa fit::TranslationalFit
   * @sa fit::PositionUtils
   * @sa utils::Rmsd
   */
  class RotationalFit{
    Reference *d_ref;

    // not implemented
    RotationalFit();
    RotationalFit(const RotationalFit &);
    RotationalFit &operator=(const RotationalFit &);

  public:
    // construct Reference for reference molecule, then
    // construct the RotationalFit.
    /**
     * RotationalFit constructor. It takes a reference to which a system 
     * can be fitted
     */
    RotationalFit(Reference *);
    /**
     * RotationalFit constructor taking an atom specifier
     */
    RotationalFit(utils::AtomSpecifier &);
    /**
     * RotationalFit deconstructor
     */
    ~RotationalFit();
    /**
     * Method to fit your System to the Reference
     */
    void fit(gcore::System *)const;
    /**
     * Method to fit your System to the Reference
     */
    void fit(utils::AtomSpecifier &, utils::AtomSpecifier &)const;

    /**
     * accessor to the reference;
     */
    Reference * getReference() { return d_ref; }
    
    struct Exception: public gromos::Exception{
      Exception(const std::string &what): 
	gromos::Exception("RotationalFit",what){}
    };
    
  };

}

#endif    
