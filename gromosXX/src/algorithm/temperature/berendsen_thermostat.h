/**
 * @file berendsen_thermostat.h
 * berendsen thermostat.
 */

#ifndef INCLUDED_TEMPERATURE_BERENDSEN_H
#define INCLUDED_TEMPERATURE_BERENDSEN_H

namespace algorithm
{
  
  /**
   * @class Berendsen_Thermostat
   * the Berendsen thermostat.
   */
  class Berendsen_Thermostat : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Berendsen_Thermostat() : Algorithm("BerendsenThermostat") {}

    /**
     * Destructor.
     */
    virtual ~Berendsen_Thermostat() {}
    
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
    
  private:

  };
  
} // algorithm

#endif
