/**
 * @file perturbed_dihedral.h
 * a perturbed dihedral.
 */

#ifndef INCLUDED_PERTURBED_DIHEDRAL_H
#define INCLUDED_PERTURBED_DIHEDRAL_H

namespace simulation
{
  /**
   * @class Perturbed_Dihedral
   * holds the perturbed dihedral information.
   */
  class Perturbed_Dihedral : public Dihedral
  {
  public:
    Perturbed_Dihedral(Dihedral &a, int B_type) 
      : Dihedral(a), B_type(B_type) {};
    int B_type;
  };
  
} // simulation

#endif
