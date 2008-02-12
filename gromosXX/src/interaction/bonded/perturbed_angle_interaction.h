/**
 * @file perturbed_angle_interaction.h
 * perturbed angle interaction.
 */

#ifndef INCLUDED_PERTURBED_ANGLE_INTERACTION
#define INCLUDED_PERTURBED_ANGLE_INTERACTION

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
   * @class Perturbed_Angle_Interaction
   * calculates the perturbed angle interactions.
   */
  class Perturbed_Angle_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Angle_Interaction(Angle_Interaction & angle_interaction)
      : Interaction("PerturbedAngle"),
	m_interaction(angle_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Angle_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      if (!quiet)
	os << "Perturbed bond angle interaction\n";
      return 0;
    };
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Angle_Interaction & m_interaction;
  };
  
} // interaction

#endif
