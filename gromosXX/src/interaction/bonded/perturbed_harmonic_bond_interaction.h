/**
 * @file perturbed_harmonic_bond_interaction.h
 * perturbed harmonic bond interaction.
 */

#ifndef INCLUDED_PERTURBED_HARMONIC_BOND_INTERACTION
#define INCLUDED_PERTURBED_HARMONIC_BOND_INTERACTION

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
   * @class Perturbed_Harmonic_Bond_Interaction
   * calculates the perturbed bond interactions (harmonic).
   */
  template<typename t_interaction_spec>
  class Perturbed_Harmonic_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Harmonic_Bond_Interaction(
       Harmonic_Bond_Interaction<t_interaction_spec> &bond_interaction)
      : Interaction("PerturbedHarmonicBond"),
	m_interaction(bond_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Harmonic_Bond_Interaction() {}

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Harmonic_Bond_Interaction<t_interaction_spec> & m_interaction;

  };
  
} // interaction

// template methods
#include "perturbed_harmonic_bond_interaction.cc"

#endif
