/**
 * @file energy_average.h
 * storage of the average energies and fluctuations.
 */

#ifndef INCLUDED_ENERGY_AVERAGE_H
#define INCLUDED_ENERGY_AVERAGE_H

namespace simulation
{
  /**
   * @class Energy_Average
   * storage of energy averages and fluctuations.
   */
  class Energy_Average
  {
  public:
    /**
     * Constructor.
     */
    Energy_Average();
    /**
     * set to zero.
     */
    void zero();
    /**
     * resize the arrays.
     */
    void resize(size_t const energy_goups, size_t const multi_baths = 0);
    /**
     * update from the calculated energies per step.
     * To allow non constant timesteps, the energies are
     * weighted by the current step size.
     */
    void update(Energy const &e, double const dt);
    /**
     * get the average energy and fluctuations
     */
    void average(Energy &energy, Energy &fluctuation);

  private:
    /**
     * the average energies.
     */
    Energy m_average;
    /**
     * the squared averages.
     */
    Energy m_square_average;
    /**
     * the time.
     */
    double m_time;

  }; // Energy_Average
  
} // simulation

#include "energy_average.tcc"

#endif

  
