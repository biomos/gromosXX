/**
 * @file periodicity.h
 * periodic boundary conditions (triclinic)
 */

#ifndef INCLUDED_PERIODICITY_H
#define INCLUDED_PERIODICITY_H

#include "math.h"

namespace math
{

  /**
   * @class Periodicity
   * the periodic boundary condition functions.
   */
  template<boundary_enum b>
  class Periodicity : public Boundary_Implementation<b>
  {
  public:
    /**
     * Constructor.
     * If b is any no specific code will be generated,
     * otherwise one can specify optimized code.
     */
    Periodicity(boundary_enum boundary = b);
    /**
     * puts a vector into the box (centered at (0, 0, 0).
     */
    void put_into_box(Vec &v)const;
    /**
     * puts a vector into the box centered at (Kx/2, Ly/2, Mz/2).
     */
    void put_into_positive_box(Vec &v)const;
  private:
  };
  
} // math

// inline functions
#include "periodicity.tcc"

#endif
  
  
    
