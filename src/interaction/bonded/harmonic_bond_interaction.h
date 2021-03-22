/**
 * @file harmonic_bond_interaction.h
 * harmonic bond interaction.
 */

#ifndef INCLUDED_HARMONIC_BOND_INTERACTION_H
#define INCLUDED_HARMONIC_BOND_INTERACTION_H

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
   * @class harmonic_bond_interaction
   * calculates the bond interactions (harmonic).
   */
  class Harmonic_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Harmonic_Bond_Interaction() : Interaction("HarmonicBond") {}
    
    /**
     * Destructor.
     */
    virtual ~Harmonic_Bond_Interaction() {}
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
      // os << "Harmonic bond interaction\n";
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
