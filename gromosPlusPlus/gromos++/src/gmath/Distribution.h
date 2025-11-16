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

// gmath_Distribution

#ifndef INCLUDED_GMATH_DISTRIBUTION
#define INCLUDED_GMATH_DISTRIBUTION

#include <iostream>
#include <vector>

#include "../gromos/Exception.h"

namespace gmath{
  
  class Vec;
  
  /**
   * Class Distribution
   * A class that calculates a distibution for a series of values
   *
   * The user has to specify an upper and lower bound and number
   * of grid points. After adding all the values a distribution can
   * be written out
   *
   * @class Distribution
   * @author B.C. Oostenbrink
   * @ingroup gmath
   * @sa gmath::stat
   */
  class Distribution{
  protected:
    double d_begin, d_end, d_step, d_sum;
    int d_nsteps, d_num;
    std::vector<int> d_count;

  public:
    /**
     * Distribution constructor
     *
     * @param begin lower bound of the distribution
     * @param end   upper bound of the distribution
     * @param nsteps number of grid points for the distribution
     */
    Distribution(double begin=0, double end=1, int nsteps=100);
    /**
     * Distribution copy constructor
     */
    Distribution(Distribution const &d);
  
    /**
     * Distribution deconstructor
     */
    ~Distribution(){}
    // Methods
    /**
     * Method to add a value to the distribution. Values outside the lower 
     * and upper bound are ignored
     * @param value The value to add
     * @return the same value is also returned
     */
    double add(const double value);
    /**
     * Method to add a vector to the distribution. Values outside the lower 
     * and upper bound are ignored
     * @param v The value to add
     * @return the same value is also returned
     */
    double add(gmath::Vec const & v)
    {
      throw Exception("vector distributions not implemented");
    }
  
    /**
     * Method to get the bin number of a given value, if that 
     * value were added to the distribution. Returns -1 if the
     * value lies outside the distribution range
     * @param value The value to hypothetically add
     * @return the bin number the value would be added in
     */
    int getbin(const double value);
    /**
     * Method to determine whether a given value lies within the distribution
     * range
     * @param value The value to hypothetically add
     * @return true or false depending on whether value lies in range
     */  
    bool inrange(const double value);
    /**
     * Method to write the complete distribution to an output stream
     * @param os an output stream (e.g. cout)
     */   
    void write(std::ostream &os)const;  
    /**
     * Method to write the complete distribution to an output stream
     * in normalized form (such that probability density integrates to 1).
     * @param os an output stream (e.g. cout)
     */
    void write_normalized(std::ostream &os)const;
    /**
     * Method to calculate the average of the values that have been added to 
     * the distribution. Values that were outside the range are not part of
     * the distribution and do not contribute to this average.
     * return The average
     */
    double ave()const;
    /**
     * Method to calculate the value x for which the distribution D is maximal,
     * i.e. D(x) = y_max. The function returns the x value, not y_max.
     * In case the distribution has to values which are maximal (D(x1) = D(x2) = y_max)
     * x1 is returned.
     */
    double maxValAt()const;
    /**
     * Method to calculate the root mean square deviation of the values that 
     * have been added. The calculation is not carried out on the values 
     * themselves, but is based on the number of elements in all the bins of 
     * the distribution.
     * @return an estimate of the rmsd
     */
    double rmsd()const;
    /**
     * Method to clear the distribution.
     */
    void clear();
  
    // Accessors
    /**
     * Accessor to obtain the number of elemenths in the i-th bin of the 
     * distribution
     * @param i the i-th bin
     * @return the number of elements in this bin
     */
    int operator[](int i)const;
    /**
     * Accessor that returns the value that corresponds to the i-th bin
     * @param i the value of the i-th bin is returned
     * @return the value that correspond to this bin (middle value)
     */
    double value(int i)const;
    /**
     * Accessor to get the number of values that have been added to the
     * distribution, values outside the lower or upper bounds are not counted.
     * @return the number of values
     */
    int nVal()const;
    /**
     * Accessor that returns the number of grid points of the distribution
     */
    int nSteps()const;
  
    // Exceptions
    struct Exception: public gromos::Exception{
      Exception(const std::string& what): 
	gromos::Exception("Distribution", what){}
    };

  };
  // inline functions & free operators
  inline double Distribution::ave()const{
    return d_sum/d_num;
  }
  inline int Distribution::operator[](int i)const{
    return d_count[i];
  }
  inline double Distribution::value(int i)const{
    return d_begin+(i+0.5)*d_step;
  }
  inline int Distribution::nVal()const
  {
    return d_num;
  }
  
  inline int Distribution::nSteps()const
  {
    return d_nsteps;
  }
  
}

#endif







