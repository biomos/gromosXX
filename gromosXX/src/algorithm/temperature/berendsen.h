/**
 * @file temperature/berendsen.h
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
  class Berendsen_Thermostat
  {
  public:
    /**
     * Constructor.
     */
    Berendsen_Thermostat();
    /**
     * apply the temperature scaling
     * or just calculate the temperature (if tau=-1).
     * @param sim the simulation.
     * @param dt the time step.
     */
    template<typename t_simulation>
    void apply(t_simulation &sim, double const dt);
    
  private:

  };
  
} // algorithm

#include "berendsen.tcc"

#endif
