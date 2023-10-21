/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
   * @param[in,out] pos the positions to transform
   * @param[in] forward true for truncoct to triclinic and false for vice versa
   */
  void truncoct_triclinic(math::VArray & pos, bool forward);
}

#endif	/* INCLUDED_TRANSFORMATION_H */
