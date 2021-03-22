/**
 * @file angle_interaction.h
 * angle interaction.
 */

#ifndef INCLUDED_ANGLE_INTERACTION_H
#define INCLUDED_ANGLE_INTERACTION_H

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
   * @class Angle_Interaction
   * calculates the angle interactions.
   */
  class Angle_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Angle_Interaction() : Interaction("Angle") {}
    /**
     * Destructor.
     */
    virtual ~Angle_Interaction() {}
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
      // os << "Bond angle (cosine) interaction\n";
      return 0;
    };
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * angle interaction parameter.
     */
    std::vector<angle_type_struct> const & parameter()const { return m_parameter; }
    /**
     * angle interaction parameter.
     */
    std::vector<angle_type_struct> & parameter() { return m_parameter; }
    
  protected:
    std::vector<angle_type_struct> m_parameter;
    
  };
  
} // interaction

#endif
