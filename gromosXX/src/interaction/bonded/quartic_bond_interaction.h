/**
 * @file quartic_bond_interaction.h
 * quartic bond interaction.
 */

#ifndef INCLUDED_QUARTIC_BOND_INTERACTION_H
#define INCLUDED_QUARTIC_BOND_INTERACTION_H

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
   * @class Quartic_bond_interaction
   * calculates the bond interactions (quartic).
   */
  class Quartic_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Quartic_Bond_Interaction() : Interaction("QuarticBond") {};
    /**
     * Destructor.
     */
    virtual ~Quartic_Bond_Interaction() {};
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology &topo,
				       configuration::Configuration &conf,
				       simulation::Simulation &sim);
    /*
     * the bond type parameter.
     */
    std::vector<bond_type_struct> const & parameter()const {return m_parameter; }
    
    /**
     * the bond type parameter.
     */
    std::vector<bond_type_struct> & parameter() { return m_parameter; }
    
  protected:
    std::vector<bond_type_struct> m_parameter;
    
  };
  
} // interaction

#endif
