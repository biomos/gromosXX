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
     * (lattice vector K)
     */
    static int const K = 0;
    /**
     * L index into the box.
     * (lattice vector L)
     */
    static int const L = 1;
    /**
     * M index into the box.
     * (lattice vector M)
     */
    static int const M = 2;

    /**
     * @struct shift_struct
     * the shift vectors.
     * used to shift any position within the box
     * to its 26 neighbours
     */
    struct shift_struct
    {
      int cell[3];
      math::Vec pos;
    };

    /**
     * Constructor.
     * @arg b the computational box
     */
    Boundary_Implementation(math::Box const & b);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     * this is the vector from v1 - v2 (gromos definition),
     * not the position of v2 closest to v1 (gromos++)
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * Get the box components of v.
     * (multipliers of K, L and M)
     */
    void box_components(Vec const &v, Vec & n)const;
    /**
     * get the computational box.
     * it's a vector of vectors K, L and M
     */
    Box const box()const;
    /**
     * get a box element.
     * returns element d2 (x, y or z) of lattice vector d1 (K, L or M)
     */
    double box(unsigned int d1, unsigned int d2)const;
    /**
     * the shifts over the periodic images.
     * shifting by lattice vectors
     * -13 .. 13
     * corresponding to (-1,-1,-1) (-1,-1,0) ... (0,0,0) ... (1,1,1)
     */
    shift_struct & shift(unsigned int i);
    /**
     * the shifts over the periodic images.
     * const accessor
     */
    shift_struct const & shift(unsigned int i)const;
    /**
     * recalculate the shift vectors.
     * necessary after box size changes (pressure coupling)
     * updates also cell index shifts
     */
    void recalc_shift_vectors(unsigned int num_cells[3]);
    /**
     * recalculate the shift vectors.
     * necessary after box size changes (pressure coupling)
     */
    void recalc_shift_vectors();

  protected:
    /**
     * reference to the box of the configuration class
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
   * the Periodicity class for the rectangular case.
   */
  template<>
  class Boundary_Implementation<math::rectangular>
  {
  public:
    /**
     * lattice vector K index
     */
    static int const K = 0;
    /**
     * lattice vector L index
     */
    static int const L = 1;
    /**
     * lattice vector M index
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
     * @arg b computational box
     */
    Boundary_Implementation(Box const & b);
    /**
     * Get the nearest image of v1 in respect to v2 (v1 - v2).
     * gromos96 convenction
     * @arg nim is delivered with the vector connecting v1 and v2
     * (not the position of v2 closest to v1, like in gromos++)
     */
    void nearest_image(Vec const &v1, Vec const &v2, Vec &nim)const;
    /**
     * Get the box components of v.
     * calculate the multipliers of K, L and M to represent vector v
     */
    void box_components(Vec const &v, Vec & n)const;
    /**
     * get the box.
     * const box accessor
     */
    Box const & box()const;
    /**
     * get a box element.
     * element d2 (x,y,z) of lattice vector d1 (K,L,M)
     */
    double box(unsigned int d1, unsigned int d2)const;
    /**
     * accessor to shift vectors over the periodic images.
     */
    shift_struct & shift(unsigned int i);
    /**
     * const accessor to shift vectors over the periodic images.
     */
    shift_struct const & shift(unsigned int i)const;
    /**
     * recalculate the shift vectors.
     * (after pressure coupling)
     */
    void recalc_shift_vectors(unsigned int num_cells[3]);

    /**
     * recalculate the shift vectors.
     * (after pressure coupling)
     */
    void recalc_shift_vectors();

  protected:
    /**
     * reference to the system::box.
     */
    Box const & m_box;
    /**
     * half box (rectangular)
     */
    math::Vec m_half_box;
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
      /**
       * shift in cell indices
       * (for grid based pairlists)
       */
      int cell[3];
      /**
       * shift vector
       */
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
     * const box accessor.
     */
    Box const & box()const;
    /**
     * get a box element.
     */
    double box(unsigned int d1, unsigned int d2)const;
    /**
     * the shifts over the periodic images.
     */
    shift_struct & shift(unsigned int i);
    /**
     * the shifts over the periodic images.
     */
    shift_struct const & shift(unsigned int i)const;
    /**
     * recalculate the shift vectors.
     */
    void recalc_shift_vectors(unsigned int num_cells[3]);
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

// template methods have to be visible
#include "boundary_implementation.cc"

#endif
