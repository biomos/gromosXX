/**
 * @file  transformation.h
 * calculate the matix transformation between coordinate systems 
 */

#ifndef INCLUDED_TRANSFORMATION_H
#define	INCLUDED_TRANSFORMATION_H

namespace math {
  /**
   * calculate the rotation matrix R from the box
   */
  math::Matrixl rmat(math::Box const & box);
    /**
   * calculate the rotation matrix R from a matrix
   */
  math::Matrix rmat(math::Matrix const & box);
   /**
   * calculate the rotation matrix R from the angles.
   */
  math::Matrixl rmat(double const & phi, double const & theta,
        double const & psi);
  /**
   * calculate the transformation matrix S.
   */
  math::Matrixl smat(math::Box const & box, math::boundary_enum const b);
  /**
   * calculate the inverse of the transformation matrix S.
   */
  math::Matrixl sinvmat(math::Box const & box, math::boundary_enum const b);
  /**
   * calculate the transformation matrix M.
   */
  math::Matrixl mmat(math::Box const & box, math::boundary_enum const b);
  /**
   * calculate the inverse of the transformation matrix M.
   */
  math::Matrixl minvmat(math::Box const & box, math::boundary_enum const b);
  /**
   * the rotation matrix for trunoct to triclinic 
   * @param[in] forward true for truncoct to triclinic and false for vice versa
   */
  math::Matrix truncoct_triclinic_rotmat(bool forward);
  /**
   * the box for trunoct to triclinic
   * @param[in] forward true for truncoct to triclinic and false for vice versa
   */
  void truncoct_triclinic_box(math::Box & box, bool forward);
  /**
   * transform a truncated octahedral box to a triclinic box and vice versa
   * @param[inout] box the box to transform
   * @param[inout] pos the positions to transform
   * @param[in] forward true for truncoct to triclinic and false for vice versa
   */
  void truncoct_triclinic(math::VArray & pos, bool forward);
}

#endif	/* INCLUDED_TRANSFORMATION_H */
