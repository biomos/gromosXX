/**
 * @file harm_angle_interaction.h
 * harmonic angle interaction.
 */

#ifndef INCLUDED_HARM_ANGLE_INTERACTION_H
#define INCLUDED_HARM_ANGLE_INTERACTION_H

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
   * @class Harm_Angle_Interaction
   * calculates the harmonic angle interactions.
   */
  class Harm_Angle_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Harm_Angle_Interaction() : Interaction("HarmAngle") {}
    /**
     * Destructor.
     */
    virtual ~Harm_Angle_Interaction() {}
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
      // os << "Harmonic bond angle interaction\n";
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
