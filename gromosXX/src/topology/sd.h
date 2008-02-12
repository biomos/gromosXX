/**
 * @file sd.h
 * stochastic dynamics variables
 */

#ifndef INCLUDED_SD_H
#define INCLUDED_SD_H

namespace topology
{
  /**
   * @struct stochastic_struct
   * holds stochastic dynamics variables
   */
  struct stochastic_struct
  {
    math::SArray gamma;
    math::SArray c1, c2, c3, c4, c5, c6, c7, c8, c9;
    
    void resize(int size)
    {
      gamma.resize(size);
      c1.resize(size);
      c2.resize(size);
      c3.resize(size);
      c4.resize(size);
      c5.resize(size);
      c6.resize(size);
      c7.resize(size);
      c8.resize(size);
      c9.resize(size);
    }
  };
}

#endif
  
