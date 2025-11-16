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

// gcore_Remd.h
#ifndef INCLUDED_REMD
#define INCLUDED_REMD

#include <cassert>
#include <vector>

#include "../gmath/Vec.h"

namespace gcore{
  /**
   * Class Remd
   * Purpose: contains REMD data
   *
   * The REMD class contains REMD data. Not a lot is done with it, but it is
   * useful if programs can read this and write it out again...
   *
   * @class Remd
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::System
   */
  class Remd{

  private:
    int d_id, d_run, d_Ti, d_li, d_Tj, d_lj, d_reeval;
    double d_temp, d_lambda;

  public:
    // constructurs
    /**
     * Remd constructor
     */
    Remd():
      d_id(0),d_run(0),d_Ti(0),d_li(0),d_Tj(0),d_lj(0),d_reeval(0),
      d_temp(0.0),d_lambda(0.0){};

    Remd(int id, int run, int Ti, int li, int Tj, int lj, int reeval,
	 double temp, double lambda):
      d_id(id),d_run(run),d_Ti(Ti),d_li(li),d_Tj(Tj),d_lj(lj),d_reeval(reeval),
      d_temp(temp),d_lambda(lambda){};

    /**
     * Remd copy constructor
     * @param r Remd to be copied
     */
    Remd(const Remd &r):
      d_id(r.d_id),d_run(r.d_run),d_Ti(r.d_Ti),d_li(r.d_li),d_Tj(r.d_Tj),
      d_lj(r.d_lj),d_reeval(r.d_reeval),d_temp(r.d_temp),d_lambda(r.d_lambda){};

    /**
     * Assignment operator
     */
    Remd &operator=(const Remd &b);
    
    // accessors
    /**
     * Accessor, returns the indentity of the remd run
     */
    int &id();
    /**
     * Accessor, returns the indentity of the remd run as a const
     */
    int id()const;
    /**
     * Accessor, returns the run-number of the remd run
     */
    int &run();
    /**
     * Accessor, returns the run-number of the remd run as a const
     */
    int run()const;
    /**
     * Accessor, returns the temperature index of state i in the remd run
     */
    int &Ti();
    /**
     * Accessor, returns the temperature index of state i in the remd run as a const
     */
    int Ti()const;
    /**
     * Accessor, returns the lambda index of state i in the remd run
     */
    int &li();   
    /**
     * Accessor, returns the lambda index of state i in the remd run as a const
     */
    int li()const;   
    /**
     * Accessor, returns the temperature index of state j in the remd run
     */
    int &Tj();
    /**
     * Accessor, returns the temperature index of state j in the remd run as a const
     */
    int Tj()const;
    /**
     * Accessor, returns the lambda index of state j in the remd run
     */
    int &lj();
    /**
     * Accessor, returns the lambda index of state j in the remd run as a const
     */
    int lj()const;
    /**
     * Accessor, returns whether it is a reeval
     */
    int &reeval();
    /**
     * Accessor, returns whether it is a reeval as a const
     */
    int reeval()const;
    /**
     * Accessor, returns the temperature of the run
     */
    double &temperature();
    /**
     * Accessor, returns the temperature of the run as a const
     */
    double temperature()const;
    /**
     * Accessor, returns the current lambda of the run
     */
    double &lambda();
    /**
     * Accessor, returns the current lambda of the run as a const
     */
    double lambda()const;
    
  };
  inline gcore::Remd &Remd::operator=(const gcore::Remd &b)
  {
    if(this!=&b){
      d_id=b.d_id;
      d_run=b.d_run;
      d_Ti=b.d_Ti;
      d_li=b.d_li;
      d_Tj=b.d_Tj;
      d_lj=b.d_lj;
      d_reeval=b.d_reeval;
      d_temp=b.d_temp;
      d_lambda=b.d_lambda;
    }
    return *this;
  }

  inline int &Remd::id(){ return d_id; }
  inline int &Remd::run(){ return d_run; }
  inline int &Remd::Ti(){ return d_Ti; }
  inline int &Remd::li(){ return d_li; }   
  inline int &Remd::Tj(){ return d_Tj; }
  inline int &Remd::lj(){ return d_lj; }
  inline int &Remd::reeval(){ return d_reeval;}
  inline double &Remd::temperature() { return d_temp; }
  inline double &Remd::lambda(){ return d_lambda; }
  inline int Remd::id()const{ return d_id; }
  inline int Remd::run()const{ return d_run; }
  inline int Remd::Ti()const{ return d_Ti; }
  inline int Remd::li()const{ return d_li; }   
  inline int Remd::Tj()const{ return d_Tj; }
  inline int Remd::lj()const{ return d_lj; }
  inline int Remd::reeval()const{ return d_reeval;}
  inline double Remd::temperature() const{ return d_temp; }
  inline double Remd::lambda()const{ return d_lambda; }
} /*namespace*/

#endif
