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
    /**
     * @struct angle_struct
     * angle information.
     */
    struct angle_struct
    {
      int i;
      int j;
      int k;
      int type;
    };

    /**
     * @class iterator
     * iterator over the angles.
     */
    class iterator
    {
    public:
      iterator(std::vector<angle_struct> &bi);
      bool eol();
      bool neol();
      void operator++();
      int i();
      int j();
      int k();
      int type();
    protected:
      std::vector<angle_struct>::const_iterator m_angle_it;
      std::vector<angle_struct>::const_iterator m_angle_end;
    };
      
    iterator begin();
    void add(int i, int j, int k, int type);

  private:      
    std::vector<angle_struct> m_angle_information;
      
  };
	  
  
} // simulation

#include "angle.tcc"

#endif
