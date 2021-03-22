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
  class Perturbed_Dihedral_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Dihedral_Interaction(Dihedral_Interaction & dihedral_interaction)
      : Interaction("PerturbedDihedral"),
	m_interaction(dihedral_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Dihedral_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
   
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Dihedral_Interaction & m_interaction;
  };
  
} // interaction

#endif
