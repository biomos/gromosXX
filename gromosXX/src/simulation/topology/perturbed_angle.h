/**
 * @file perturbed_angle.h
 * a perturbed angle.
 */

#ifndef INCLUDED_PERTURBED_ANGLE_H
#define INCLUDED_PERTURBED_ANGLE_H

namespace simulation
{
  /**
   * @class Perturbed_Angle
   * holds the perturbed bond information.
   */
  class Perturbed_Angle : public Angle
  {
  public:
    Perturbed_Angle(Angle &a, int B_type) : Angle(a), B_type(B_type) {};
    int B_type;
  };
  
} // simulation

#endif
