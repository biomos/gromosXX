/**
 * @file perturbed_improper_dihedral.h
 * a perturbed improper dihedral.
 */

#ifndef INCLUDED_PERTURBED_IMPROPER_DIHEDRAL_H
#define INCLUDED_PERTURBED_IMPROPER_DIHEDRAL_H

namespace simulation
{
  /**
   * @class Perturbed_Angle
   * holds the perturbed bond information.
   */
  class Perturbed_Improper_Dihedral : public Improper_Dihedral
  {
  public:
    Perturbed_Improper_Dihedral(Improper_Dihedral &a, int B_type) 
      : Improper_Dihedral(a), B_type(B_type) {};
    int B_type;
  };
  
} // simulation

#endif
