/**
 * @file gmath.h
 * mathematical definitions.
 */

#ifndef INCLUDED_MATH_H
#define INCLUDED_MATH_H

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

#include <random/normal.h>

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <typeinfo>

#ifdef COMPILER_GCC
#include <cxxabi.h>
#endif
/**
 * @namespace math
 * namespace that contains mathematical functions
 * using Blitz++ (www.oonumerics.org/blitz)
 */
namespace math
{

  using namespace blitz;
  BZ_USING_NAMESPACE(blitz)

  using namespace ranlib;

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
  typedef blitz::TinyMatrix<double, 3U, 3U> Matrix;

  /**
   * Box.
   */
  typedef blitz::TinyVector< blitz::TinyVector<double, 3>, 3> Box;

  /**
   * @enum boundary_enum
   * boundary condition
   */
  enum boundary_enum{
    /**
     * vacuum.
     */
    vacuum = 0,
    /**
     * triclinic box
     */
    triclinic = 1,
    /**
     * non-specialized version
     */
    any = 2
  };
  
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

  /**
   * h bar.
   */
  extern double h_bar;
  
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

/**
 * blitz dot product
 */
template<typename T, int N>
inline T dot(TinyVector<T,N> const &v1, TinyVector<T,N> const &v2)
{
  return blitz::dot(v1, v2);
}

BZ_DECLARE_FUNCTION2_RET(dot, double);

/**
 * blitz cross product
 */
inline math::Vec cross(math::Vec const &v1, math::Vec const & v2)
{
  return blitz::cross(v1, v2);
}

BZ_DECLARE_FUNCTION2_RET(cross, math::Vec);

/**
 * blitz abs2
 */
template<typename T, int N>
inline T abs2(TinyVector<T,N> v)
{
  return sum(sqr(v));
}

BZ_DECLARE_FUNCTION_RET(abs2, double);

} // math

#endif

