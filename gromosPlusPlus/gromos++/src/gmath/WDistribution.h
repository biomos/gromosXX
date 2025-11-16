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

// gmath_WDistribution

#ifndef INCLUDED_GMATH_WDISTRIBUTION
#define INCLUDED_GMATH_WDISTRIBUTION

#include <iostream>
#include <vector>

#include "Stat.h"
#include "../gromos/Exception.h"

namespace gmath{
  
  class Vec;
  
  /**
   * Class WDistribution
   * A class that calculates a distibution for a series of values with
   * different weights.
   *
   * The user has to specify an upper and lower bound and number
   * of grid points. After adding all the values a distribution can
   * be written out
   *
   * @class WDistribution
   * @author @ref cc 
   * @ingroup gmath
   * @sa gmath::stat
   */
  class WDistribution : public Distribution {
    // all weights are now floating point numbers and not unity like in Distribution
    Stat<double> d_num_d_stat;
    double d_num_d;
    std::vector<Stat<double> > d_count_d_stat;
    std::vector<double> d_count_d;

  public:
    /**
     * WDistribution constructor
     *
     * @param begin lower bound of the distribution
     * @param end   upper bound of the distribution
     * @param nsteps number of grid points for the distribution
     */
    WDistribution(double begin=0, double end=1, int nsteps=100);
    /**
     * WDistribution copy constructor
     */
    WDistribution(WDistribution const &d);

    /**
     * WDistribution deconstructor
     */
    ~WDistribution(){}
    /**
     * Method to add a value with its corresponding weight to the distribution.
     * Values outside the lower
     * and upper bound are ignored
     * @param value The value to add
     * @param weight The logarithm of the corresponding unnormalized weight
     * @return the same value is also returned
     */
    double add(const double value, const double weight);
    /**
     * Method to write the complete distribution to an output stream
     * @param os an output stream (e.g. cout)
     */
    void write(std::ostream &os);
    /**
     * Method to write the complete distribution to an output stream
     * in normalized form.
     * @param os an output stream (e.g. cout)
     */
    void write_normalized(std::ostream &os);
    /**
     * Method to calculate the average of the values that have been added to
     * the distribution. Values that were outside the range are not part of
     * the distribution and do not contribute to this average.
     * return The average
     */
    double ave()const{
      throw WDistribution::Exception("ave() not implemented.");
    }
    /**
     * Method to calculate the root mean square deviation of the values that
     * have been added. The calculation is not carried out on the values
     * themselves, but is based on the number of elements in all the bins of
     * the distribution.
     * @return an estimate of the rmsd
     */
    double rmsd()const{
      throw WDistribution::Exception("rmsd() not implemented.");
    }
    /**
     * Method to clear the distribution.
     */
    void clear(){
      throw WDistribution::Exception("clear() not implemented.");
    }

    // Accessors
    /**
     * Accessor to obtain the number of elemenths in the i-th bin of the
     * distribution
     * @param i the i-th bin
     * @return the number of elements in this bin
     */
    int operator[](int i)const{
      throw WDistribution::Exception("operator[] not implemented.");
    }

     // Exceptions
    struct Exception: public gromos::Exception{
      Exception(const std::string& what):
	gromos::Exception("WDistribution", what){}
    };
  };
  
}

#endif







