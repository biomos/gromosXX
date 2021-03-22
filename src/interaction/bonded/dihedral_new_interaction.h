/**
 * @file dihedral_new_interaction.h
 * dihedral interaction.
 */

#ifndef INCLUDED_DIHEDRAL_NEW_INTERACTION_H
#define INCLUDED_DIHEDRAL_NEW_INTERACTION_H

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
   * @class Dihedral_new_Interaction
   * calculates the dihedral interactions.
   */
  class Dihedral_new_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Dihedral_new_Interaction() : Interaction("Dihedral") {}
    /**
     * Destructor.
     */
    virtual ~Dihedral_new_Interaction() {}

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
      // os << "Dihedral interaction\n";
      return 0;
    };
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * the angle type parameters.
     */
    std::vector<dihedral_type_struct> const & parameter()const { return m_parameter; }
    /**
     * the angle type parameters.
     */
    std::vector<dihedral_type_struct> & parameter() { return m_parameter; }
    
  protected:
    std::vector<dihedral_type_struct> m_parameter;

    /**
     * calculate nearest minimum
     */
    // double _calculate_nearest_minimum(double phi, int m, double pd);
    
  };
  
} // interaction

#endif
