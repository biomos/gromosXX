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
   * implements periodic boundary conditions.
   */
  template<math::boundary_enum b>
  class Boundary_Implementation;

  /**
   * @class Boundary_Implementation<vacuum>
   * Specialized version for vacuum.
   */
  template<>
  class Boundary_Implementation<vacuum>
  {
  public:
    /**
     * K index into the box.
     */
    static int const K = 0;
    /**
     * L index into the box.
     */
    static int const L = 1;
    /**
     * M index into the box.
     */
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
     */
    Boundary_Implementation(math::Box const & b);
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
    /**
     * get the box.
     */
    Box const box()const;
    /**
     * get a box element.
     */
    double const box(size_t const d1, size_t const d2)const;

    /**
     * the shifts over the periodic images.
     */
    shift_struct & shift(size_t const i);

    /**
     * the shifts over the periodic images.
     */
    shift_struct const & shift(size_t const i)const;

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors(size_t const num_cells[3]);

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors();

  protected:
    /**
     * reference to the system::box.
     */
    Box const & m_box;

    /**
     * the shift vectors.
     */
    shift_struct m_shift[27];

  };

  /*
   * @class Boundary_Implementation
   * implements the specific functions of
   * the Periodicity class.
   */
  template<>
  class Boundary_Implementation<math::rectangular>
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
     */
    Boundary_Implementation(Box const & b);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * Get the box components of v.
     */
    void box_components(Vec const &v, Vec & n)const;
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
     * the shifts over the periodic images.
     */
    shift_struct & shift(size_t const i);

    /**
     * the shifts over the periodic images.
     */
    shift_struct const & shift(size_t const i)const;

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors(size_t const num_cells[3]);

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors();

  protected:
    /**
     * reference to the system::box.
     */
    Box const & m_box;
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
     */
    Boundary_Implementation(Box const & b);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * Get the box components of v.
     */
    void box_components(Vec const &v, Vec & n)const;

    /**
     * get the box.
     */
    Box const & box()const;
    /**
     * get a box element.
     */
    double const box(size_t const d1, size_t const d2)const;

    /**
     * the shifts over the periodic images.
     */
    shift_struct & shift(size_t const i);

    /**
     * the shifts over the periodic images.
     */
    shift_struct const & shift(size_t const i)const;

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors(size_t const num_cells[3]);

    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors();

  protected:
    /**
     * reference to the system::box.
     */
    Box const & m_box;
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

#include "boundary_implementation.cc"

#endif
