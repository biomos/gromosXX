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
 * @file TranslationalFit.h
 */

#ifndef INCLUDED_FIT_TRANSLATIONALFIT
#define INCLUDED_FIT_TRANSLATIONALFIT


namespace gcore{
  class System;
}

namespace gmath{
  class Vec;
}

namespace fit{
  class TranslationalFit_i;
  class Reference;

  enum centre_enum { cog, com };

  /**
   * Class TranslationalFit
   * A class that performs a translational fit of one system with 
   * respect to a reference
   *
   * @class TranslationalFit
   * @author R. Buergi
   * @ingroup fit
   * @sa RotationalFit
   * @sa Reference
   */
  class TranslationalFit{
    TranslationalFit_i *d_this;

    // not implemented
    TranslationalFit();
    TranslationalFit(const TranslationalFit &);
    TranslationalFit &operator=(const TranslationalFit &);

  public:
    // construct Reference for reference molecule, then
    // construct the TranslationalFit.
    /**
     * TranslationalFit constructor.
     * @arg Reference reference to which your system will be fitted
     * @arg centre cog or com : centre of geometry or centre of mass.
     *      default is centre of geometry
     * 
     */
    TranslationalFit(Reference *ref, centre_enum centre = fit::cog);
    /**
     * TranslationalFit deconstructor
     */
    ~TranslationalFit();
    
    //    void fitToCom(gcore::System *)const;
    /**
     * Method that fits the System to the Reference
     */
    void fit(gcore::System *)const;
    /**
     * Accessor that returns the centre of mass of the system including the 
     * weight defined in the Reference
     */
    const gmath::Vec &com()const;
    /**
     * Accessor that returns the centre of geometry of the system, including
     * the weights defined in the Reference
     */
    const gmath::Vec &cog()const;
    
  };

}

#endif    
