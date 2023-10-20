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
 * @file bs_vector.h
 * 
 * The special vectors used for B&S-LEUS, where the dimension is arbitrary long
 * and periodicities are taken into account.
 */

#ifndef BS_VECTOR_H
#define	BS_VECTOR_H

namespace util {
    
  /**
   * @class BS_Vector
   * 
   * Implementation for an n-dimensional vector holding the coordinates
   * of the subspace
   */
  class BS_Vector : public std::vector<double> {
  public:
    /**
     * return the length squared
     * @return squared length
     */
    double abs2();
    /**
     * this - subtrahend = result
     * @param[in]    subtrahend
     * @param[in,out] result
     */
    void minus(const BS_Vector &subtrahend, BS_Vector &result);
    /**
     * scale the vector by scalar
     * @param[in] scalar
     */
    void scale(const double scalar);
    /**
     * Normalize the vector to length 1
     * @return the original length
     */
    double normalize();
    /**
     * Set every entry to zero.
     */
    void nullify();
    /**
     * Multiply the vector by scalar
     * @param[in] scalar
     * @return the scaled vector
     */
    BS_Vector operator*(const double scalar);
    /**
     * Add two vectors together
     * @param[in] summand
     * @return the sum
     */
    BS_Vector operator +(const BS_Vector &summand);
    /**
     * Subtract two vectors from each other
     * @param[in] subtrahend
     * @return the difference
     */
    BS_Vector operator -(const BS_Vector &subtrahend);
    /**
     * Add summand to a vector
     * @param summand
     */
    void operator +=(const BS_Vector &summand);
    /**
     * The dot product with other:
     *      dot(self, other)
     * @param other vector
     * @return the dot product
     */
    double dot(const BS_Vector &other);
    /**
     * Create a BS_Vector from a vector of doubles
     * @param values
     */
    void create(std::vector<double> &values);
    BS_Vector operator=(std::vector<double> &values){create(values); return *this;}
   /**
     * Creates an output for the Vector.
     * @return the output
     */
    std::string str();
  };
}


#endif	/* BS_VECTOR_H */

