/**
 * @file berendsen_thermostat.h
 * berendsen thermostat.
 */

#ifndef INCLUDED_TEMPERATURE_BERENDSEN_H
#define INCLUDED_TEMPERATURE_BERENDSEN_H

#include "thermostat.h"

namespace algorithm
{
  
  /**
   * @class Berendsen_Thermostat
   * the Berendsen thermostat.
   */
  class Berendsen_Thermostat : public Thermostat
  {
  public:
    /**
     * Constructor.
     */
    Berendsen_Thermostat() : Thermostat("BerendsenThermostat") {}

    /**
     * Destructor.
     */
    virtual ~Berendsen_Thermostat() {}
    
    /**
     * initialise
     */
    virtual int init
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     std::ostream & os = std::cout,
     bool quiet = false
     );

    /**
     * apply the temperature scaling
     * for baths with tau=-1 nothing is done.
     * the kinetic energy can not be calculated here, because
     * later on SHAKE might be applied.
     * @param topo the Topology
     * @param conf the Configuration
     * @param sim the Simulation
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    /**
     * calculate the scaling factors.
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param immediate if true rescales the velocities to immediately satisfy
     * the given reference temperature (strong coupling).
     */
    void calc_scaling(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      bool immediate = false);

  private:

  };
  
} // algorithm

#endif
