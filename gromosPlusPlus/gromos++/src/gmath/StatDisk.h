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

// gmath_STATDISK

#ifndef INCLUDED_GMATH_STATDISK
#define INCLUDED_GMATH_STATDISK

#include <vector>
#include <fstream>
#include <string>

namespace gmath
{
  /**
   * Class StatDisk
   * A class to perform some basic statistics on a series of numbers
   *
   * This class allows one to store a series of numbers and calculate 
   * the average, rmsd and an error estimate
   *
   * The data is stored in a scratch file and not in the memory. This 
   * (and lack of random access and distributions) is the difference to
   * Stat
   *
   * @class StatDisk
   * @author N. Schmid
   * @ingroup gmath
   * @sa Stat
   */
  template<typename T>
  class StatDisk
  {
    mutable std::vector<int> d_blocksize;
    std::vector<T> d_vals;
    int d_counter;
    mutable T d_ave,d_msd, d_ee;
    mutable bool d_avedone, d_msddone, d_eedone;
    std::string d_file;
    std::ofstream d_out;
    size_t d_autoflush;
      
  public:
    /**
     * StatDisk constructor
     * @param file the scratch file
     * @param autoflush the maximum size of the internal buffer
     */
    StatDisk(std::string file, size_t autoflush = 1000);
    /**
     * StatDisk constructor
     * @param autoflush the maximum size of the internal buffer
     */
    StatDisk(size_t autoflush = 1000);
    /**
     * StatDisk destructor
     */
    ~StatDisk();
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
    T msd();
    /**
     * Method to calculate the square root of the mean square deviation
     * @return root-mean-square-deviation
     * @sa msd
     */
    T rmsd();
    /**
     * Method to calculate the average over the series. Internally, we
     * determine whether a new calculation is required or if we can just 
     * return the previously calculated value.
     * @return The average
     */
    T ave();
    /**
     * Method to calculate the average over only part of the series.
     * @param b first index of the series
     * @param e last index of the series. The average is calculated 
     *          for(i=b; i<e; i++)
     * @return The average of this range
     */
    T subave(int b, int e);
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
    T ee();
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
    T min();
    /**
     * Accessor that returns the maximum value of the values stored so
     * far
     * requires that operator> is defined for type T
     */
    T max();
    /**
     * flush the buffer to the file
     */
    void flush();
    /**
     * open file file
     */
    void open(std::string file);
  };
  
  template<typename T>
  inline int StatDisk<T>::n()const
  {
    return d_counter;
  }
  
}

#include "StatDisk.cc"
#endif
