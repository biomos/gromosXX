/**
 * @file perturbed_harmonic_bond_interaction.h
 * perturbed harmonic bond interaction.
 */

#ifndef INCLUDED_PERTURBED_HARMONIC_BOND_INTERACTION
#define INCLUDED_PERTURBED_HARMONIC_BOND_INTERACTION

namespace interaction
{
  /**
   * @class Perturbed_Harmonic_Bond_Interaction
   * calculates the perturbed bond interactions (harmonic).
   */
  template<typename t_simulation>
  class Perturbed_Harmonic_Bond_Interaction : public Interaction<t_simulation>
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Harmonic_Bond_Interaction(
       harmonic_bond_interaction<t_simulation> &bond_interaction);
    /**
     * Destructor.
     */
    virtual ~Perturbed_Harmonic_Bond_Interaction();
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual void calculate_interactions(t_simulation &sim);
    
  protected:
    harmonic_bond_interaction<t_simulation> &m_bond_interaction;
  };
  
} // interaction

// template methods
#include "perturbed_harmonic_bond_interaction.tcc"

#endif
