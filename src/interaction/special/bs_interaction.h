/** 
 * @file bs_interaction.h
 * 
 * Contains the B&S-LEUS Algorithm Interface.
 */

#ifndef BS_INTERACTION_H
#define	BS_INTERACTION_H

namespace interaction {
  /**
   * @class BS_Interaction Implements the B&S-LEUS algorithm of Halvor Hansen
   */
  class BS_Interaction : public Interaction {
  public:

    /**
     * Constructor.
     */
    BS_Interaction() : Interaction("B&S-LEUS") {}
    /**
     * Destructor.
     */
    virtual ~BS_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os = std::cout,
            bool quiet = false);
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);
  protected:
  };
}

#endif	/* BS_INTERACTION_H */

