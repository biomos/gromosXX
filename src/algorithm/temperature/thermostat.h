/**
 * @file thermostat.h
 * thermostat.
 */

#ifndef INCLUDED_TEMPERATURE_THERMOSTAT_H
#define INCLUDED_TEMPERATURE_THERMOSTAT_H

namespace algorithm
{
  
  /**
   * @class Thermostat
   * thermostat abstract base
   */
  class Thermostat : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Thermostat(std::string s) : Algorithm(s) {}

    /**
     * Destructor.
     */
    virtual ~Thermostat() {}
    
    /**
     * rescale the velocities.
     */
    void scale(topology::Topology & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation & sim);
    
  private:
    
  };
  
} // algorithm

#endif
