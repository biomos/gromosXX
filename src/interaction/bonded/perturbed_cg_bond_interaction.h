/**
 * @file perturbed_cg_bond_interaction.h
 * perturbed cg bond interaction.
 */

#ifndef INCLUDED_PERTURBED_CG_BOND_INTERACTION
#define INCLUDED_PERTURBED_CG_BOND_INTERACTION

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
   * @class Perturbed_CG_Bond_Interaction
   * calculates the perturbed bond interactions (cg).
   */
  class Perturbed_DP_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_DP_Bond_Interaction(
       DP_Bond_Interaction &bond_interaction)
      : Interaction("PerturbedDPBond"),
	m_interaction(bond_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_DP_Bond_Interaction() {}

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
      // os << "Perturbed cg bond interaction\n";
      return 0;
    };

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    DP_Bond_Interaction & m_interaction;

  };
  
} // interaction

#endif
