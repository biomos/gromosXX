/**
 * @file gmath.h
 * mathematical definitions.
 */

#ifndef INCLUDED_MATH_H
#define INCLUDED_MATH_H

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

#include <vector>
#include <map>
#include <set>

/**
 * @namespace math
 * namespace that contains mathematical functions
 * using Blitz++ (www.oonumerics.org/blitz)
 */
namespace math
{

  BZ_USING_NAMESPACE(blitz)

  /**
   * 3 dimensional vector.
   */
  typedef blitz::TinyVector<double, 3U> Vec;
  /**
   * Array of 3D vectors.
   */
  typedef blitz::Array<Vec, 1>         VArray;
  /**
   * Array of scalars.
   */
  typedef blitz::Array<double, 1>      SArray;
  /**
   * Matrix.
   */
  typedef blitz::TinyVector< blitz::TinyVector<double, 3>, 3> Matrix;
  
  /**
   * a small number.
   */
  const double epsilon = 0.000000000001;

  /**
   * Pi
   */
  const double Pi = 3.1415926535897932384626433;

  /**
   * Boltzmann constant.
   */
  const double k_Boltzmann = 0.00831441;
  
#ifndef NDEBUG
  /**
   * module debug level.
   */
  extern int debug_level;
#endif


/**
 * provide comparision operators for the blitz TinyVector.
 * they should be implemented by blitz, but i cannot get
 * them to work?!
 */
inline bool operator==(math::Vec &t1, math::Vec &t2)
{
  bool b = true;
  for(int i=0; i<3; ++i)
    if (t1(i) != t2(i)) b = false;
  return b;
}

/**
 * != operator
 */
inline bool operator!=(math::Vec &t1, math::Vec &t2)
{
  return !(t1 == t2);
}

} // math

#include "boundary_implementation.h"
#include "periodicity.h"

#endif

