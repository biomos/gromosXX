/**
 * @file boundary_checks.h
 * checks for boxes and boundary conditions
 */

#ifndef INCLUDED_BOUNDARY_CHECKS_H
#define	INCLUDED_BOUNDARY_CHECKS_H

namespace math {
  /**
   * checks whether the box is large enough for a given cutoff
   * @param box the box
   * @param b the boundary conditions
   * @param cutoff the cutoff
   */
  bool boundary_check_cutoff(math::Box const & box, math::boundary_enum const b,
          double cutoff);
}

#endif	/* _BOUNDARY_CHECKS_H */

