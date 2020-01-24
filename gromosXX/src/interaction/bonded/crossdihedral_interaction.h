/**
 * @file crossdihedral_interaction.h
 * crossdihedral interaction.
 */

#ifndef INCLUDED_CROSSDIHEDRAL_INTERACTION_H
#define INCLUDED_CROSSDIHEDRAL_INTERACTION_H

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
   * @class Crossdihedral_Interaction
   * calculates the crossdihedral interactions.
   */
  class Crossdihedral_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Crossdihedral_Interaction() : Interaction("Crossdihedral") {}
    /**
     * Destructor.
     */
    virtual ~Crossdihedral_Interaction() {}

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
      // os << "Crossdihedral interaction\n";
      return 0;
    };
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:

    /**
     * calculate nearest minimum
     */
    // double _calculate_nearest_minimum(double phi, int m, double pd);
    
  };
  
} // interaction

#endif
