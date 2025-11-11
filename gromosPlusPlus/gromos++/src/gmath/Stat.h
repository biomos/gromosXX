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

// gmath_STAT

#ifndef INCLUDED_GMATH_STAT
#define INCLUDED_GMATH_STAT

#include <vector>

#include "Distribution.h"
#include "../gromos/Exception.h"

namespace gmath
{
  /**
   * Class Stat
   * A class to perform some basic statistics on a series of numbers
   *
   * This class allows one to store a series of numbers and calculate 
   * the average, rmsd and an error estimate
   *
   * @class Stat
   * @author B.C. Oostenbrink
   * @ingroup gmath
   */
  template<typename T>
  class Stat
  {
    mutable std::vector<int> d_blocksize;
    std::vector<T> d_vals;
    int d_counter;
    mutable T d_ave, d_lnexpave, d_msd, d_ee;
    mutable bool d_avedone, d_lnexpavedone, d_msddone, d_eedone, d_distdone;
    gmath::Distribution d_dist;
      
  public:
    /**
     * Stat constructor
     */
    Stat();
    /**
     * Stat destructor
     */
    ~Stat(){}
    /**
     * Method to add another value to the series
     * @param val the value to add
     */
    void addval(T val);
    /**
     * Method to calculate (or return) the mean square deviation 
     * of the series. Within the class, we keep track of whether 
     * anything has changed since the previous calculation to determine
     * if a new calculation is needed
     * @return mean-square-deviation
     */
    T msd()const;
    /**
     * Method to calculate the square root of the mean square deviation
     * @return root-mean-square-deviation
     * @sa msd
     */
    T rmsd()const;
    /**
     * Method to calculate the average over the series. Internally, we
     * determine whether a new calculation is required or if we can just 
     * return the previously calculated value.
     * @return The average
     */
    T ave()const;
    /**
     * Method to calculate the exponential average ln<exp(X)> over
     * the series. Internally, we determine whether a new calculation is
     * required or if we can just return the previously calculated value.
     * @return The average
     */
    T lnexpave()const;
    /**
     * Method to calculate the average over only part of the series.
     * @param b first index of the series
     * @param e last index of the series. The average is calculated 
     *          for(i=b; i<e; i++)
     * @return The average of this range
     */
    T subave(int b, int e)const;
    /**
     * Method to calculate the exponential average ln<exp(X)> over only part of the series.
     * @param b first index of the series
     * @param e last index of the series. The average is calculated
     *          for(i=b; i<e; i++)
     * @return The average of this range
     */
    T lnexpsubave(int b, int e)const;
    /**
     * Method to calculate the average ln|<Xexp(Y)>| over the series.
     * @return ln|<Xexp(Y)>|
     */
    static T lnXexpave(Stat<T> X, Stat<T> Y, int &sign);
    /**
     * Method to calculate the covariance Cov(X,Y)=\<XY\>-\<X\>\<Y\>
     * @param X a Stat object
     * @param Y another Stat object
     * @return  the covariance
     */
    static T covariance(Stat<T> X, Stat<T> Y);
    /**
     * Method to calculate ln{Cov(exp(X),exp(Y))}=ln{<exp(X)exp(Y)>-<exp(X)><exp(Y)>}
     * @param X a Stat object
     * @param Y another Stat object
     * @return the logarithm of the covariance and the sign of it.
     */
    static T lnexpcovariance(Stat<T> X, Stat<T> Y, int & sign);
    /**
     * Method to calculate an error estimate for the series.
     *
     * The error estimation is based on a method described by Alan and 
     * Tildesley. The series is devided into blocks, for which the average
     * is calculated. The rmsd of these averages is then calculated for 
     * different block sizes. An extrapolation to infinite block size then 
     * gives the error estimate.
     * @return The error estimate
     */
    T ee()const;
    /**
     * Method to calculate the statistical inefficiency as described in
     * J. Chem. Theory Comput. 2007, 3, 26-41. The implementation closely
     * follows the Python code provided by the authors.
     */
    static T stat_ineff(Stat<T> X, Stat<T> Y);
    /**
     * Method to calculate the logarithm of the statistical inefficiency as described in
     * J. Chem. Theory Comput. 2007, 3, 26-41. The implementation closely
     * follows the Python code provided by the authors. It returns the logarithm of the
     * statistical inefficiency of x and y with X=ln(x) and Y=ln(y).
     */
    static T lnexp_stat_ineff(Stat<T> X, Stat<T> Y);
    /**
     * Accessor that returns the stored values
     * @param i the i-th value is returned
     * @return the value that was stored
     */
    T val(int i);
    /**
     * Accessor to return the number of elements that have been stored 
     * so far
     * @return the number of values that are stored in the class
     */
    int n()const;
    /**
     * Accessor that returns the minimum value of the values stored so 
     * far
     * requires that operator< is defined for type T
     */
    T min()const;
    /**
     * Accessor that returns the maximum value of the values stored so
     * far
     * requires that operator> is defined for type T
     */
    T max()const;
    /**
     * Accessor that returns a pointer to the vector containing the data
     */
    std::vector<T> const & data() const;
    /**
     * Accessor that returns a pointer to a Distribution of the data
     *
     * The stat class also contains a Distribution of the data. 
     * Call dist_init to define the bounds of the Distribution.
     * @return a pointer to a gmath::Distribution
     */
    gmath::Distribution const & distribution()const;
    /**
     * Initializes a Distribution
     *
     * The function dist_init rescales the Distribution to the specified
     * bounds and puts all values so far in the distribution.
     * Values that are added later, are also added to the distribution.
     * @param lower lower bound of the distribution
     * @param upper upper bound of the distribution
     * @param nsteps number of grid points for the distribution
     * @param periodic assumes (upper - lower) to be a full period and maps all values periodically into this interval
     */
    gmath::Distribution const & dist_init(double lower, double upper, int nsteps, bool periodic = false);
    /**
     * Initializes a Distribution
     *
     * The function dist_init rescales the Distribution to use the current
     * minima and maxima as bounds.
     * @param nsteps number of grid points for the distribution
     */
    gmath::Distribution const & dist_init(int nsteps);
    /**
     * Subtract the average from every data point
     */
    void subtract_average();
      
  };

  //inline functions
  template<typename T>
  inline T Stat<T>::val(int i)
  {
    return d_vals[i];
  }
  
  template<typename T>
  inline int Stat<T>::n()const
  {
    return d_counter;
  }

  template<typename T>
  inline gmath::Distribution const & Stat<T>::dist_init(int nsteps)
  {
    if (this->max() == this->min()){
      const double dd = this->max() * 0.01;
      return dist_init(this->min() - dd, this->max() + dd, nsteps);
    }
    return this->dist_init(this->min(), this->max()+1E-10, nsteps);
  }
  
}

#include "Stat.cc"

#endif
