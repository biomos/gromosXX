/**
 * @file external_interaction.h
 * external interaction
 */

#ifndef INCLUDED_EXTERNAL_INTERACTION_H
#define INCLUDED_EXTERNAL_INTERACTION_H

namespace interaction
{
  /**
   * @class External_Interaction
   * implements an external interaction
   */
  class External_Interaction : public interaction::Interaction
  {
  public:
    /**
     * Constructor.
     */
    External_Interaction() : Interaction("External"), cg_topo(NULL), cg_conf(NULL) {}
    
    /**
     * Destructor.
     */
    virtual ~External_Interaction() {}
    
    /**
     * External interaction:
     * - virtual graining interaction
     */
    virtual int calculate_interactions(topology::Topology &topo, 
				       configuration::Configuration &conf,
				       simulation::Simulation &sim);
    /**
     * initialisation
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream & os = std::cout,
		     bool quiet = false)
    {
      os << "External interaction\n";
      return 0;
    }
    

    void set_coarsegraining(topology::Topology & cg_topo,
			    configuration::Configuration & cg_conf,
			    simulation::Simulation & cg_sim)
    {
      this->cg_topo = &cg_topo;
      this->cg_conf = &cg_conf;
      this->cg_sim = & cg_sim;
    }
    
  private:
    /**
     * coarse grained topology
     * (virtual atom (grains) definition)
     */
    topology::Topology * cg_topo;
    /**
     * coarse grained configuration
     * (forces on virtual atoms (grains))
     */
    configuration::Configuration * cg_conf;
    /**
     * coarse grained simulation
     * (input parameters might be different!)
     */
    simulation::Simulation * cg_sim;
    
  };
  
} // interaction

#endif

