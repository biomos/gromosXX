/**
 * @file perturbed_bond.h
 * a perturbed bond.
 */

#ifndef INCLUDED_PERTURBED_BOND_H
#define INCLUDED_PERTURBED_BOND_H

namespace simulation
{
  /**
   * @class Perturbed_Bond
   * holds the perturbed bond information.
   */
  class Perturbed_Bond : public Bond
  {
  public:
    Perturbed_Bond(Bond &b, int B_type) : Bond(b), B_type(B_type) {};
    int B_type;
  };
  
} // simulation

#endif
