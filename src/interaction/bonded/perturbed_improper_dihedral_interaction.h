/**
 * @file perturbed_improper_dihedral_interaction.h
 * perturbed improper dihedral interaction.
 */

#ifndef INCLUDED_PERTURBED_IMPROPER_DIHEDRAL_INTERACTION
#define INCLUDED_PERTURBED_IMPROPER_DIHEDRAL_INTERACTION

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
   * @class Perturbed_Improper_Dihedral_Interaction
   * calculates the perturbed angle interactions.
   */
  class Perturbed_Improper_Dihedral_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Improper_Dihedral_Interaction
      (Improper_Dihedral_Interaction & improper_dihedral_interaction)
	: Interaction("PerturbedImproperDihedral"),
	  m_interaction(improper_dihedral_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Improper_Dihedral_Interaction() {}

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
      // os << "Perturbed improper dihedral interaction\n";
      return 0;
    };
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Improper_Dihedral_Interaction & m_interaction;

  };
  
} // interaction

#endif
