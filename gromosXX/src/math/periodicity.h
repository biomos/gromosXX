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
    triclinic = 1
  };

  /**
   * @class periodicity
   * the periodic boundary condition functions.
   */
  template<boundary_enum b>
  class periodicity;
  
  /**
   * @class periodicityvacuum>
   * Specialized version for vacuum.
   */
  template<>
  class periodicity<vacuum>
  {
  public:
    periodicity(Matrix &box);
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim);
  private:
    Matrix &m_box;
  };
  
  /**
   * @class peridicity<triclinic>
   * specialized version for triclinic boundary conditions.
   */
  template<>
  class periodicity<triclinic>
  {
  public:
    periodicity(Matrix &box);
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim);
  private:
    Matrix &m_box;
  };
  
} // math

// inline functions
#include "periodicity.tcc"

#endif
  
  
    
