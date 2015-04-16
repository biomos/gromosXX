/**
 * @file cg_bond_interaction.h
 * cg bond interaction.
 */

#ifndef INCLUDED_CG_BOND_INTERACTION_H
#define INCLUDED_CG_BOND_INTERACTION_H

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
   * @class cg_bond_interaction
   * calculates the bond interactions (cg).
   */
  class DP_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    DP_Bond_Interaction() : Interaction("DPBond") {}
    
    /**
     * Destructor.
     */
    virtual ~DP_Bond_Interaction() {}
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
      // os << "CG bond interaction\n";
      return 0;
    };

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * the bond parameter.
     */
    std::vector<bond_type_struct> const & parameter()const { return m_parameter;}
    /**
     * the bond parameter.
     */
    std::vector<bond_type_struct> & parameter() { return m_parameter; }
    
  protected:
    std::vector<bond_type_struct> m_parameter;
    
  };
  
} // interaction

#endif
