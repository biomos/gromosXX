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

#ifndef NDEBUG
  /**
   * module debug level.
   */
  extern int debug_level;
#endif

} // math

#include "periodicity.h"

#endif

