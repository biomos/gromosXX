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
     * const perturbed bonds.
     */
    std::vector<Perturbed_Bond> const & bonds()const;

    /**
     * perturbed bonds.
     */
    std::vector<Perturbed_Bond> & bonds();

    /**
     * const perturbed angles.
     */
    std::vector<Perturbed_Angle> const & angles()const;
    
    /**
     * perturbed angles.
     */
    std::vector<Perturbed_Angle> & angles();
    
    /**
     * perturbed atoms accessor.
     */
    std::map<size_t, Perturbed_Atom> & atoms();

    /**
     * const perturbed atoms accessor.
     */
    std::map<size_t, Perturbed_Atom> const & atoms()const;
    
  private:
    /**
     * the perturbed bonds.
     */
    std::vector<Perturbed_Bond> m_bond;

    /**
     * the perturbed angles.
     */
    std::vector<Perturbed_Angle> m_angle;
    
    /**
     * the perturbed atoms.
     */
    std::map<size_t, Perturbed_Atom> m_atom;
    
  };
  
} // simulation

#include "perturbed_solute.tcc"

#endif
