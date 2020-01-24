/**
 * @file nemd.h
 * nemd
 */

#ifndef INCLUDED_NEMD_INTERACTION_H
#define INCLUDED_NEMD_INTERACTION_H

namespace interaction
{
  /**
   * @class nemd_interaction
   * non-equilibrium molecular dynamics
   */
  class NEMD_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    NEMD_Interaction() : Interaction("NEMD") {}

    /**
     * Destructor.
     */
    virtual ~NEMD_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false)
    {
      os << "NEMD interaction\n";
      return 0;
    }
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

  protected:

  };

} // interaction

#endif

