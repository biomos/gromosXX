/**
 * @file gmath.h
 * mathematical definitions.
 */

#ifndef INCLUDED_MATH_H
#define INCLUDED_MATH_H

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

#include <vector>

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
  typedef blitz::TinyVector<double, 3> Vec;
  /**
   * Array of 3D vectors.
   */
  typedef blitz::Array<Vec, 1>         VArray;
  /**
   * Array of scalars.
   */
  typedef blitz::Array<double, 1>      SArray;
  
  
} // math

#endif

