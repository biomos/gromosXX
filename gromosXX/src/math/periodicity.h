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
    triclinic = 1,
    /**
     * non-specialized version
     */
    any = 2
  };

  /**
   * @class Periodicity
   * the periodic boundary condition functions.
   */
  template<boundary_enum b>
  class Periodicity
  {
  public:
    Periodicity(Matrix &box, boundary_enum boundary = b);
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    void boundary(boundary_enum const b);
    boundary_enum const boundary()const;
  private:
    Matrix &m_box;
    boundary_enum m_boundary;
  };
  
  /**
   * @class Periodicity<vacuum>
   * Specialized version for vacuum.
   */
  template<>
  class Periodicity<vacuum>
  {
  public:
    Periodicity(Matrix &box, boundary_enum boundary = vacuum);
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    void boundary(boundary_enum const b);
    boundary_enum const boundary()const;
  private:
    Matrix &m_box;
  };
  
  /**
   * @class Periodicity<triclinic>
   * specialized version for triclinic boundary conditions.
   */
  template<>
  class Periodicity<triclinic>
  {
  public:
    Periodicity(Matrix &box, boundary_enum boundary = triclinic);
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    void boundary(boundary_enum const b);
    boundary_enum const boundary()const;
  private:
    Matrix &m_box;
  };
  
} // math

// inline functions
#include "periodicity.tcc"

#endif
  
  
    
