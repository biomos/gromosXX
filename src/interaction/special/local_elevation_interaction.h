/**
 * @file local_elevation_interaction.h
 * local elevation umbrella samping interaction
 */

#ifndef INCLUDED_LOCAL_ELEVATION_INTERACTION_H
#define	INCLUDED_LOCAL_ELEVATION_INTERACTION_H

namespace interaction {

  /**
   * @class xray_restraint_interaction
   * calculates the xray restraining interaction
   */ class Local_Elevation_Interaction : public Interaction {
  public:

    /**
     * Constructor.
     */
    Local_Elevation_Interaction() : Interaction("Local Elevation") {}
    /**
     * Destructor.
     */
    virtual ~Local_Elevation_Interaction() {}

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
} // interaction
#endif	/* INCLUDED_LOCAL_ELEVATION_INTERACTION_H */

