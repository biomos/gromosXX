/**
 * @file improper_dihedral.h
 * the improper dihedral topology class.
 */

#ifndef INCLUDED_IMPROPER_DIHEDRAL_H
#define INCLUDED_IMPROPER_DIHEDRAL_H

namespace simulation
{
  /**
   * @class Improper_Dihedral
   * holds improper dihedral information.
   */
  class Improper_Dihedral
  {
  public:
    Improper_Dihedral(int i, int j, int k, int l, int type)
      : i(i), j(j), k(k), l(l), type(type) {};

    int i;
    int j;
    int k;
    int l;
    int type;
  };
	  
  
} // simulation

#endif
