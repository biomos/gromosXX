/**
 * @file boundary_implementation.h
 * periodic boundary conditions (triclinic)
 * nearest image implementation.
 */

#ifndef INCLUDED_BOUNDARY_IMPLEMENTATION_H
#define INCLUDED_BOUNDARY_IMPLEMENTATION_H

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
   * @class Boundary_Implementation
   * implements the specific functions of
   * the Periodicity class.
   */
  template<boundary_enum b>
  class Boundary_Implementation
  {
  public:
    /**
     * Constructor.
     * @param box refernce to the box member of system.
     * @param boundary is the boundary condition.
     */
    Boundary_Implementation(Matrix &box, boundary_enum boundary = b);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * set the boundary condition.
     */
    void boundary(boundary_enum const b);
    /**
     * get the boundary condition.
     */
    boundary_enum const boundary()const;
  protected:
    /**
     * reference to the system::box.
     */
    Matrix &m_box;
    /**
     * the boundary condition.
     */
    boundary_enum m_boundary;
  };
  
  /**
   * @class Boundary_Implementation<vacuum>
   * Specialized version for vacuum.
   */
  template<>
  class Boundary_Implementation<vacuum>
  {
  public:
    /**
     * Constructor.
     * @param box refernce to the box member of system.
     * @param boundary is the boundary condition.
     */
    Boundary_Implementation(Matrix &box, boundary_enum boundary = vacuum);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * set the boundary condition.
     */
    void boundary(boundary_enum const b);
    /**
     * get the boundary condition.
     */
    boundary_enum const boundary()const;
  protected:
    /**
     * reference to the system::box.
     */
    Matrix &m_box;
  };
  
  /**
   * @class Boundary_Implementation<triclinic>
   * specialized version for triclinic boundary conditions.
   */
  template<>
  class Boundary_Implementation<triclinic>
  {
  public:
    /**
     * Constructor.
     * @param box refernce to the box member of system.
     * @param boundary is the boundary condition.
     */
    Boundary_Implementation(Matrix &box, boundary_enum boundary = triclinic);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * set the boundary condition.
     */
    void boundary(boundary_enum const b);
    /**
     * get the boundary condition.
     */
    boundary_enum const boundary()const;
  protected:
    /**
     * reference to the system::box.
     */
    Matrix &m_box;
  };
  
}

#include "boundary_implementation.tcc"

#endif
