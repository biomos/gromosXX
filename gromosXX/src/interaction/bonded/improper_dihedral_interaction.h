/**
 * @file improper_dihedral_interaction.h
 * improper dihedral interaction.
 */

#ifndef INCLUDED_IMPROPER_DIHEDRAL_INTERACTION_H
#define INCLUDED_IMPROPER_DIHEDRAL_INTERACTION_H

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
   * @class Improper_Dihedral_Interaction
   * calculates the improper dihedral interactions.
   */
  class Improper_Dihedral_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Improper_Dihedral_Interaction() : Interaction("ImproperDihedral") {}
    /**
     * Destructor.
     */
    virtual ~Improper_Dihedral_Interaction() {}
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
      // os << "Improper dihedral interaction\n";
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
    std::vector<improper_dihedral_type_struct> const & parameter()const { return m_parameter; }

    /**
     * the angle type parameters.
     */
    std::vector<improper_dihedral_type_struct> & parameter() { return m_parameter; }
    
  protected:
    std::vector<improper_dihedral_type_struct> m_parameter;
    
  };
  
} // interaction

#endif
