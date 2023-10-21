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
 * @file sd.h
 * stochastic dynamics variables
 */

#ifndef INCLUDED_SD_H
#define INCLUDED_SD_H

namespace topology
{
  /**
   * @struct stochastic_struct
   * holds stochastic dynamics variables
   */
  struct stochastic_struct
  {
    math::SArray gamma;
    math::SArray c1, c2, c3, c4, c5, c6, c7, c8, c9;
    /**
     * resize the variables
     */
    void resize(int size)
    {
      gamma.resize(size);
      c1.resize(size);
      c2.resize(size);
      c3.resize(size);
      c4.resize(size);
      c5.resize(size);
      c6.resize(size);
      c7.resize(size);
      c8.resize(size);
      c9.resize(size);
    }
    
    /**
     * clear the variables
     */
    void clear() {
      gamma.clear();
      c1.clear();
      c2.clear();
      c3.clear();
      c4.clear();
      c5.clear();
      c6.clear();
      c7.clear();
      c8.clear();
      c9.clear();      
    }
    
    /**
     * constructor
     */
    stochastic_struct() :
    gamma(0), c1(0), c2(0), c3(0), c4(0), c5(0), c6(0), c7(0), c8(0), c9(0) {}
    
    /**
     * copy constructor
     */
    stochastic_struct(const stochastic_struct& copy) {
      gamma = copy.gamma;
      c1 = copy.c1;
      c2 = copy.c2;
      c3 = copy.c3;
      c4 = copy.c4;
      c5 = copy.c5;
      c6 = copy.c6;
      c7 = copy.c7;
      c8 = copy.c8;
      c9 = copy.c9;
    }
  };
}

#endif
  
