/**
 * @file perturbed_quartic_bond_interaction.h
 * perturbed quartic bond interaction.
 */

#ifndef INCLUDED_PERTURBED_QUARTIC_BOND_INTERACTION
#define INCLUDED_PERTURBED_QUARTIC_BOND_INTERACTION

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
   * @class Perturbed_Quartic_Bond_Interaction
   * calculates the perturbed bond interactions (quartic).
   */
  class Perturbed_Quartic_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Quartic_Bond_Interaction(Quartic_Bond_Interaction &bond_interaction)
      : Interaction("PerturbedQuarticBond"),
	m_interaction(bond_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Quartic_Bond_Interaction() {}

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
      // os << "Perturbed quartic bond interaction\n";
      return 0;
    };

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Quartic_Bond_Interaction & m_interaction;

  };
  
} // interaction

#endif
