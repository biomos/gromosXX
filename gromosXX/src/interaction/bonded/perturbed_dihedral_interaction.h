/**
 * @file perturbed_dihedral_interaction.h
 * perturbed  dihedral interaction.
 */

#ifndef INCLUDED_PERTURBED_DIHEDRAL_INTERACTION
#define INCLUDED_PERTURBED_DIHEDRAL_INTERACTION

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

namespace interaction
{
  /**
   * @class Perturbed_Dihedral_Interaction
   * calculates the perturbed dihedral interactions.
   */
  template<typename t_interaction_spec>
  class Perturbed_Dihedral_Interaction 
    : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Dihedral_Interaction
      (Dihedral_Interaction<t_interaction_spec> 
       & dihedral_interaction)
	: Interaction("PerturbedDihedral"),
	  m_interaction(dihedral_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Dihedral_Interaction() {}

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Dihedral_Interaction<t_interaction_spec> & m_interaction;
  };
  
} // interaction

// template methods
#include "perturbed_dihedral_interaction.cc"

#endif
