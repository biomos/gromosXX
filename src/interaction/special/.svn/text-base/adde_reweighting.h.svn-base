/**
 * @file adde_reweighting.h
 * reweighting for adiabatic decoupling
 */

#ifndef INCLUDED_ADDE_REWEIGHTING_H
#define INCLUDED_ADDE_REWEIGHTING_H

namespace interaction
{
  /**
   * @class adde_reweighting
   */
  class Adde_Reweighting : public Interaction
  {
  public:
    /**
     * Constructor.
     */
     Adde_Reweighting(): Interaction("AddeReweighting") {}
    
    /**
     * Destructor.
     */
    virtual ~Adde_Reweighting() {}

     /**
     * init
     */
    virtual int init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os = std::cout,
            bool quiet = false) {
      conf.special().adde.evhl = 0;
      return 0;
    }
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    
  };
  
} // interaction

#endif

