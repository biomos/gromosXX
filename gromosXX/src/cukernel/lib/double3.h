/**
 * @file double3.h
 * 3D vector operations
 */

#ifndef _DOUBLE3_H
#define	_DOUBLE3_H

#ifndef HOSTDEVICE
#error "Don't include double3.h without defining HOSTDEVICE"
#else

#include "types.h"
#include <cstdlib>
/**
 * produce a double3
 * this is not required in CUDA 3.2 anymore and thus commented out...
 *
__host__ HOSTDEVICE double3 make_double3(const double & a, const double & b, const double & c){
  double3 resultat;
  resultat.x = a;
  resultat.y = b;
  resultat.z = c;
  return resultat;
}
*/

/**
 * calculates the scalar (dot) product of two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return the dot product
 */
HOSTDEVICE double dot(const double3 & a, const double3 & b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * calculates the vector (cross) product of two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return the dot product
 */
HOSTDEVICE double3 cross(const double3 & a, const double3 & b) {
  return make_double3(a.y*b.z - a.z*b.y,
                     a.z*b.x - a.x*b.z,
                     a.x*b.y - a.y*b.x);
}

/**
 * calculates the squared length of a vector
 * @param a the vector
 * @return squared length
 */
HOSTDEVICE double abs2(const double3 & a) {
  return dot(a,a);
}

/**
 * calculates the length of a vector
 * @param a the vector
 * @return length
 */
HOSTDEVICE double abs(const double3 & a) {
  return sqrt(abs2(a));
}

/**
 * operator to subtract two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return a-b
 */
HOSTDEVICE double3 operator-(const double3 & a, const double3 & b) {
  return make_double3(a.x-b.x, a.y-b.y, a.z-b.z);
}

/**
 * operator to add two vectors
 * @param[in] a first vector
 * @param[out] b second vector
 * @return a+b
 */
HOSTDEVICE double3 operator+(const double3 & a, const double3 & b) {
  return make_double3(a.x+b.x, a.y+b.y, a.z+b.z);
}

/**
 * negates a vector
 * @param[in] a the vector
 * @return -a
 */
HOSTDEVICE double3 operator-(const double3 & a) {
  return make_double3(-a.x,-a.y,-a.z);
}

/**
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
HOSTDEVICE double3 operator*(const double3 & a, double b) {
  return make_double3(a.x*b, a.y*b, a.z*b);
}

/**
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
HOSTDEVICE double3 operator*(double b, const double3 & a) {
  return make_double3(a.x*b, a.y*b, a.z*b);
}

/**
 * scales a vector (division by a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a/b
 */
HOSTDEVICE double3 operator/(const double3 & a, double b) {
  b = 1.0f / b;
  return a*b;
}

#endif

#endif	/* _DOUBLE3_H */

