/**
 * @file boundary_checks.cc
 * checks for boxes and boundary conditions (implementation)
 */

#include <stdheader.h>
#include "boundary_checks.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

bool math::boundary_check_cutoff(math::Box const & box, math::boundary_enum const b,
        double cutoff) {
  
  const double cutoff2 = 2.0 * cutoff;
  DEBUG(3, "Checking cutoff: " << cutoff);
  DEBUG(3, "Checking box: " << std::endl << math::v2s(box(0)) << std::endl << math::v2s(box(1)) << std::endl << math::v2s(box(2)));
  
  switch (b) {
    case math::vacuum :
      break;
    case math::rectangular :
    {
      if (abs(box(0)) <= cutoff2 ||
          abs(box(1)) <= cutoff2 ||
          abs(box(2)) <= cutoff2) {
        // box too small
        return false;
      }

      break;
    }
    case math::triclinic :
    {
      double a, b, c, alpha, beta, gamma, triclinicvolume;
      a = math::abs(box(0));
      b = math::abs(box(1));
      c = math::abs(box(2));

      alpha = acos(dot(box(1), box(2)) / (abs(box(1)) * abs(box(2))));
      beta = acos(dot(box(0), box(2)) / (abs(box(0)) * abs(box(2))));
      gamma = acos(dot(box(0), box(1)) / (abs(box(0)) * abs(box(1))));

      triclinicvolume = a * b * c *
              (1.0 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) - cos(gamma) * cos(gamma)
              + 2.0 * cos(alpha) * cos(beta) * cos(gamma));

      if (triclinicvolume / (a * b * sin(gamma)) <= cutoff2 ||
          triclinicvolume / (a * c * sin(beta)) <= cutoff2 ||
          triclinicvolume / (b * c * sin(alpha)) <= cutoff2) {
        return false; // too small
      }

      break;
    }
    case math::truncoct :
    {
      if (0.5 * sqrt(3.0) * abs(box(0)) <= cutoff2) {
        return false; // too small
      }
      break;
    }
    default:
      io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
  }
  
  // the box is fine
  return true;
}

