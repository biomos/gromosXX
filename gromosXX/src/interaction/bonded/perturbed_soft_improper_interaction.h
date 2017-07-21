/**
 * @file perturbed_soft_improper_interaction.h
 * perturbed soft improper dihedral interaction.
 */

#ifndef INCLUDED_PERTURBED_SOFT_IMPROPER_INTERACTION
#define INCLUDED_PERTURBED_SOFT_IMPROPER_INTERACTION

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
   * @class Perturbed_Soft_Improper_Interaction
   * calculates the perturbed soft improper dihedral interactions.
   */
  class Perturbed_Soft_Improper_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Soft_Improper_Interaction()
      : Interaction("PerturbedSoftImproper"){}
        
    /**
     * Destructor.
     */
    virtual ~Perturbed_Soft_Improper_Interaction() {}

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
    

  };
  
} // interaction

#endif
