/**
 * @file dihedral.h
 * the dihedral topology class.
 */

#ifndef INCLUDED_DIHEDRAL_H
#define INCLUDED_DIHEDRAL_H

namespace simulation
{
  /**
   * @class Dihedral
   * holds dihedral information.
   */
  class Dihedral
  {
  public:
    Dihedral(int i, int j, int k, int l, int type)
      : i(i), j(j), k(k), l(l), type(type) {};
    
    int i;
    int j;
    int k;
    int l;
    int type;
  };
	  
  
} // simulation

#endif
