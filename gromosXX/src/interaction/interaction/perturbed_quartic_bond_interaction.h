/**
 * @file perturbed_quartic_bond_interaction.h
 * perturbed quartic bond interaction.
 */

#ifndef INCLUDED_PERTURBED_QUARTIC_BOND_INTERACTION
#define INCLUDED_PERTURBED_QUARTIC_BOND_INTERACTION

namespace interaction
{
  /**
   * @class Perturbed_Quartic_Bond_Interaction
   * calculates the perturbed bond interactions (quartic).
   */
  template<typename t_simulation>
  class Perturbed_Quartic_Bond_Interaction : public Interaction<t_simulation>
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Quartic_Bond_Interaction(Quartic_bond_interaction<t_simulation> 
				       &bond_interaction);
    /**
     * Destructor.
     */
    virtual ~Perturbed_Quartic_Bond_Interaction();
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual void calculate_interactions(t_simulation &sim);
    
  protected:
    Quartic_bond_interaction<t_simulation> &m_bond_interaction;
  };
  
} // interaction

// template methods
#include "perturbed_quartic_bond_interaction.tcc"

#endif
