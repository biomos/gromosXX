/**
 * @file boundary_checks.cc
 * checks for boxes and boundary conditions (implementation)
 */

#include "../stdheader.h"
#include "boundary_checks.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

bool math::boundary_check_cutoff(math::Box const & box, math::boundary_enum const b,
        double cutoff) {
  
  const double cutoff2 = 2.0 * cutoff;
  DEBUG(12, "Checking cutoff: " << cutoff);
  DEBUG(12, "Checking box: " << std::endl << math::v2s(box(0)) << std::endl << math::v2s(box(1)) << std::endl << math::v2s(box(2)));
  
  switch (b) {
    case math::vacuum :
      DEBUG(13, "PBC: vacuum");
      break;
    case math::rectangular :
    {
      DEBUG(13, "PBC: rectangular");
      if (abs(box(0)) <= cutoff2 ||
          abs(box(1)) <= cutoff2 ||
          abs(box(2)) <= cutoff2) {
        // box too small
        DEBUG(13, "box is too small");
        return false;
      }

      break;
    }
    case math::truncoct :
    case math::triclinic :
    {
      DEBUG(13, "PBC: truncoct/triclinic");
      double a = 0.0, b = 0.0, c = 0.0, alpha = 0.0, beta = 0.0, gamma = 0.0, triclinicvolume = 0.0;
      a = math::abs(box(0));
      b = math::abs(box(1));
      c = math::abs(box(2));
      DEBUG(13, "a: " << a << " b: " << b << " c: " << c);

      alpha = acos(dot(box(1), box(2)) / (abs(box(1)) * abs(box(2))));
      beta = acos(dot(box(0), box(2)) / (abs(box(0)) * abs(box(2))));
      gamma = acos(dot(box(0), box(1)) / (abs(box(0)) * abs(box(1))));
      DEBUG(13, "al: " << alpha * 180.0 / math::Pi
              << " be: " << beta * 180.0 / math::Pi
              << " ga: " << gamma * 180.0 / math::Pi);

      triclinicvolume = a * b * c *
              sqrt(1.0 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) - cos(gamma) * cos(gamma)
              + 2.0 * cos(alpha) * cos(beta) * cos(gamma));

      DEBUG(13, "vol: " << triclinicvolume);

      const double dim1 = triclinicvolume / (a * b * sin(gamma));
      const double dim2 = triclinicvolume / (a * c * sin(beta));
      const double dim3 = triclinicvolume / (b * c * sin(alpha));
      DEBUG(13, "dims: " << dim1 << " " << dim2 << " " << dim3);

      if (dim1 <= cutoff2 || dim2 <= cutoff2 || dim3 <= cutoff2) {
        DEBUG(13, "box is too small");
        return false; // too small
      }

      break;
    }
    default:
      io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
  }
  
  // the box is fine
  DEBUG(12, "box is fine");
  return true;
}

