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
  template<typename t_interaction_spec>
  class Perturbed_Improper_Dihedral_Interaction 
    : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Improper_Dihedral_Interaction
      (Improper_Dihedral_Interaction<t_interaction_spec> 
       & improper_dihedral_interaction)
	: Interaction("PerturbedImproperDihedral"),
	  m_interaction(improper_dihedral_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Improper_Dihedral_Interaction() {}

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Improper_Dihedral_Interaction<t_interaction_spec> &m_interaction;

  };
  
} // interaction

// template methods
#include "perturbed_improper_dihedral_interaction.tcc"

#endif
