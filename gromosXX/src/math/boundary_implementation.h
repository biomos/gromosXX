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
   * @class Boundary_Implementation
   * implements the specific functions of
   * the Periodicity class.
   */
  template<boundary_enum b>
  class Boundary_Implementation
  {
  public:
    static int const K = 0;
    static int const L = 1;
    static int const M = 2;

    /**
     * @struct shift_struct
     * the shift vectors.
     */
    struct shift_struct
    {
      int cell[3];
      math::Vec pos;
    };

    /**
     * Constructor.
     * @param boundary is the boundary condition.
     */
    Boundary_Implementation(boundary_enum boundary = b);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * Get the box components of v.
     */
    void box_components(Vec const &v, Vec & n)const;
    /**
     * set the boundary condition.
     */
    void boundary_condition(boundary_enum const b);
    /**
     * get the boundary condition.
     */
    boundary_enum const boundary_condition()const;
    /**
     * get the box.
     */
    Box const & box()const;
    /**
     * get a box vector.
     */
    // Vec const & box(size_t const d)const;
    /**
     * get a box element.
     */
    double const box(size_t const d1, size_t const d2)const;

    /**
     * set the box.
     */
    void box(Box const &m);
    // not implemented
    /**
     * set the box.
     */
    // void box(Vec v1, Vec v2, Vec v3);
    
    /**
     * the volume
     */
    double volume()const;

    /**
     * the shifts over the periodic images.
     */
    shift_struct & shift(size_t const i);

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors(size_t const num_cells[3]);

  protected:
    /**
     * reference to the system::box.
     */
    Box m_box;
    /**
     * the boundary condition.
     */
    boundary_enum m_boundary;
    /**
     * the box volume.
     */
    double m_volume;
    /**
     * triclinic nearest image:
     * -(L*M) / vol
     * -(K*M) / vol
     * -(K*L) / vol
     */
    Box m_cross_K_L_M;

    /**
     * the shift vectors.
     */
    shift_struct m_shift[27];

  };
  
  /**
   * @class Boundary_Implementation<vacuum>
   * Specialized version for vacuum.
   */
  template<>
  class Boundary_Implementation<vacuum>
  {
  public:
    static int const K = 0;
    static int const L = 1;
    static int const M = 2;

    /**
     * @struct shift_struct
     * the shift vectors.
     */
    struct shift_struct
    {
      int cell[3];
      math::Vec pos;
    };

    /**
     * Constructor.
     * @param boundary is the boundary condition.
     */
    Boundary_Implementation(boundary_enum boundary = vacuum);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * Get the box components of v.
     */
    void box_components(Vec const &v, Vec & n)const;
    /**
     * set the boundary condition.
     */
    void boundary_condition(boundary_enum const b);
    /**
     * get the boundary condition.
     */
    boundary_enum const boundary_condition()const;
    /**
     * get the box.
     */
    Box const box()const;
    /**
     * get a box vector.
     */
    // Vec const box(size_t const d)const;
    /**
     * get a box element.
     */
    double const box(size_t const d1, size_t const d2)const;

    /**
     * set the box.
     */
    void box(Box const &m);
    // not implemented
    /**
     * set the box.
     */
    // void box(Vec v1, Vec v2, Vec v3);

    /**
     * the volume
     */
    double volume()const;

    /**
     * the shifts over the periodic images.
     */
    shift_struct & shift(size_t const i);

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors(size_t const num_cells[3]);

  protected:
    /**
     * reference to the system::box.
     */
    Box m_box;
    /**
     * the box volume.
     */
    double m_volume;

    /**
     * the shift vectors.
     */
    shift_struct m_shift[27];

  };
  
  /**
   * @class Boundary_Implementation<triclinic>
   * specialized version for triclinic boundary conditions.
   */
  template<>
  class Boundary_Implementation<triclinic>
  {
  public:
    static int const K = 0;
    static int const L = 1;
    static int const M = 2;

    /**
     * @struct shift_struct
     * the shift vectors.
     */
    struct shift_struct
    {
      int cell[3];
      math::Vec pos;
    };

    /**
     * Constructor.
     * @param boundary is the boundary condition.
     */
    Boundary_Implementation(boundary_enum boundary = triclinic);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * Get the box components of v.
     */
    void box_components(Vec const &v, Vec & n)const;
    /**
     * set the boundary condition.
     */
    void boundary_condition(boundary_enum const b);
    /**
     * get the boundary condition.
     */
    boundary_enum const boundary_condition()const;

    /**
     * get the box.
     */
    Box const & box()const;
    /**
     * get a box vector.
     */
    // Vec const & box(size_t const d)const;
    /**
     * get a box element.
     */
    double const box(size_t const d1, size_t const d2)const;

    /**
     * set the box.
     */
    void box(Box const &m);
    // not implemented
    /**
     * set the box.
     */
    // void box(Vec v1, Vec v2, Vec v3);

    /**
     * the volume.
     */
    double volume()const;

    /**
     * the shifts over the periodic images.
     */
    shift_struct & shift(size_t const i);

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors(size_t const num_cells[3]);

  protected:
    /**
     * reference to the system::box.
     */
    Box m_box;
    /**
     * the box volume.
     */
    double m_volume;
    /**
     * triclinic nearest image:
     * -(L*M) / vol
     * -(K*M) / vol
     * -(K*L) / vol
     */
    Box m_cross_K_L_M;

    /**
     * the shift vectors.
     */
    shift_struct m_shift[27];

  };
  
}

#include "boundary_implementation.tcc"

#endif
