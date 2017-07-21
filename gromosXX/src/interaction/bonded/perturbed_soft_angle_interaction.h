/**
 * @file perturbed_soft_angle_interaction.h
 * perturbed soft angle interaction.
 */

#ifndef INCLUDED_PERTURBED_SOFT_ANGLE_INTERACTION
#define INCLUDED_PERTURBED_SOFT_ANGLE_INTERACTION

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

namespace io{
	class IFP;
}

namespace interaction
{
  /**
   * @class Perturbed_Angle_Interaction
   * calculates the perturbed angle interactions.
   */
  class Perturbed_Soft_Angle_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Soft_Angle_Interaction(io::IFP &it);
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Soft_Angle_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) ;
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    std::vector<angle_type_struct> m_parameter;
  };
  
} // interaction

#endif
