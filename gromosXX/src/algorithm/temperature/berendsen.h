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
     * if tau=-1 nothing is done.
     * the kinetic energy can not be calculated here, because
     * later on SHAKE might be applied.
     * @param sim the simulation.
     * @param dt the time step.
     */
    template<typename t_simulation>
    void apply(t_simulation &sim, double const dt);

    /**
     * calculate the kinetic energy.
     * @param sim the simulation.
     */
    template<typename t_simulation>
    void calculate_kinetic_energy(t_simulation &sim);

    /**
     * calculate the kinetic energy lambda derivative
     * @param sim the simulation.
     */
    template<typename t_simulation>
    void calculate_kinetic_energy_lambda_derivative(t_simulation &sim);
    
  private:

  };
  
} // algorithm

#include "berendsen.tcc"

#endif
