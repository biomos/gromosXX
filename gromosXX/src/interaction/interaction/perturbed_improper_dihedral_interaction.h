/**
 * @file perturbed_improper_dihedral_interaction.h
 * perturbed improper dihedral interaction.
 */

#ifndef INCLUDED_PERTURBED_IMPROPER_DIHEDRAL_INTERACTION
#define INCLUDED_PERTURBED_IMPROPER_DIHEDRAL_INTERACTION

namespace interaction
{
  /**
   * @class Perturbed_Improper_Dihedral_Interaction
   * calculates the perturbed angle interactions.
   */
  template<typename t_simulation, typename t_interaction_spec>
  class Perturbed_Improper_Dihedral_Interaction 
    : public Interaction<t_simulation, t_interaction_spec>
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Improper_Dihedral_Interaction
      (Improper_dihedral_interaction<t_simulation, t_interaction_spec> 
       & improper_dihedral_interaction);
    /**
     * Destructor.
     */
    virtual ~Perturbed_Improper_Dihedral_Interaction();
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual void calculate_interactions(t_simulation &sim);
    
  protected:
    Improper_dihedral_interaction<t_simulation, t_interaction_spec>
         &m_improper_dihedral_interaction;
  };
  
} // interaction

// template methods
#include "perturbed_improper_dihedral_interaction.tcc"

#endif
