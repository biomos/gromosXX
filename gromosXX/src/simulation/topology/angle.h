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

    bool operator==(Angle const &a)
    {
      return (i==a.i && j==a.j && k==a.k && type==a.type);
    };
    
  };
	  
  
} // simulation

#endif
