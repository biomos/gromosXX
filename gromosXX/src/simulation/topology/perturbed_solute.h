/**
 * @file perturbed_solute.h
 * the perturbed part of the solute.
 */

#ifndef INCLUDED_PERTURBED_SOLUTE_H
#define INCLUDED_PERTURBED_SOLUTE_H

namespace simulation
{
  /**
   * @class Perturbed_Solute
   * holds the information about perturbation of
   * the solute.
   */
  class Perturbed_Solute
  {
  public:
    /**
     * perturbed bonds.
     */
    std::vector<Perturbed_Bond> & bonds();

  private:
    /**
     * the perturbed bonds.
     */
    std::vector<Perturbed_Bond> m_bond;
  };
  
} // simulation

#include "perturbed_solute.tcc"

#endif
