/**
 * @file bond.h
 * the bond topology class.
 */

#ifndef INCLUDED_BOND_H
#define INCLUDED_BOND_H

namespace simulation
{
  /**
   * @class Bond
   * holds bond information.
   */
  class Bond
  {
  public:
    Bond(int i, int j, int type) : i(i), j(j), type(type) {};
    int i;
    int j;
    int type;

    bool operator==(Bond const &b)
    {
      return (i==b.i && j==b.j && type == b.type);
    };
    
  };
	  
  
} // simulation

#endif
