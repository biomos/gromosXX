/**
 * @file perturbed_dihedral_interaction.h
 * perturbed  dihedral interaction.
 */

#ifndef INCLUDED_PERTURBED_DIHEDRAL_INTERACTION
#define INCLUDED_PERTURBED_DIHEDRAL_INTERACTION

namespace interaction
{
  /**
   * @class Perturbed_Dihedral_Interaction
   * calculates the perturbed dihedral interactions.
   */
  template<typename t_simulation>
  class Perturbed_Dihedral_Interaction 
    : public Interaction<t_simulation>
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Dihedral_Interaction
      (Dihedral_interaction<t_simulation> 
       & dihedral_interaction);
    /**
     * Destructor.
     */
    virtual ~Perturbed_Dihedral_Interaction();
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual void calculate_interactions(t_simulation &sim);
    
  protected:
    Dihedral_interaction<t_simulation>
         &m_dihedral_interaction;
  };
  
} // interaction

// template methods
#include "perturbed_dihedral_interaction.tcc"

#endif
