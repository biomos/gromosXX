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

    bool operator==(Improper_Dihedral const &id)
    {
      return (i==id.i && j==id.j && k==id.k && l==id.l && type==id.type);
    };

  };
	  
  
} // simulation

#endif
