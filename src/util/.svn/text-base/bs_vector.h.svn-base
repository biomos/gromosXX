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
     * @param[inout] result
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
     * @param[in] summand
     * @return the sum
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
     * Create a BS_Vector from a <double> vector     
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

