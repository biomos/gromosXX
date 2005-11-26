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
  class Perturbed_Harmonic_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Harmonic_Bond_Interaction(
       Harmonic_Bond_Interaction &bond_interaction)
      : Interaction("PerturbedHarmonicBond"),
	m_interaction(bond_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Harmonic_Bond_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      // if (!quiet)
      // os << "Perturbed harmonic bond interaction\n";
      return 0;
    };

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Harmonic_Bond_Interaction & m_interaction;

  };
  
} // interaction

#endif
