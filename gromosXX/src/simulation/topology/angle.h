/**
 * @file angle.h
 * the angle topology class.
 */

#ifndef INCLUDED_ANGLE_H
#define INCLUDED_ANGLE_H

namespace simulation
{
  /**
   * @class Angle
   * holds angle information.
   */
  class Angle
  {
  public:
    Angle(int i, int j, int k, int type) : i(i), j(j), k(k), type(type) {};
    
    int i;
    int j;
    int k;
    int type;
  };
	  
  
} // simulation

#endif
