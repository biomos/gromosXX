/**
 * @file volume.cc
 */

#include <util/stdheader.h>
#include "volume.h"

double math::volume(math::Box & box, math::boundary_enum b)
{
  switch (b){
    case math::vacuum:
      return 0;
    case math::rectangular:
      return box(0)(0) * box(1)(1) * box(2)(2);
    case math::triclinic:
      return math::dot(math::cross(box(0), box(1)), box(2));
    default:
      return 0;
  }
  return 0;
}
