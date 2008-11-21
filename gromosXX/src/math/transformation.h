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
}

#endif	/* INCLUDED_TRANSFORMATION_H */
