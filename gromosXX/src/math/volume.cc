/**
 * @file volume.cc
 */

#include <stdheader.h>
#include "volume.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

double math::volume(math::Box const & box, math::boundary_enum const b)
{
  switch (b){
    case math::vacuum:
      DEBUG(9, "vacuum: volume = 0");
      return 0;
    case math::rectangular:
      DEBUG(9, "rectangular: volume = "
	    <<  box(0)(0) * box(1)(1) * box(2)(2));
      return box(0)(0) * box(1)(1) * box(2)(2);
    case math::triclinic:
      DEBUG(9, "triclinic: volume = "
	    << math::dot(math::cross(box(0), box(1)), box(2)));
      return math::dot(math::cross(box(0), box(1)), box(2));
    default:
      DEBUG(9, "volume error....");
      return 0;
  }
  return 0;
}
