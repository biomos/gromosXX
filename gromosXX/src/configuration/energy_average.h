/**
 * @file energy_average.h
 * storage of the average energies and fluctuations.
 */

#ifndef INCLUDED_ENERGY_AVERAGE_H
#define INCLUDED_ENERGY_AVERAGE_H

namespace configuration
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
    void update(configuration::Energy const &e, double const dt);

    /**
     * average the pressure.
     */
    void update(math::Matrix const &pressure, double const dt);

    /**
     * get the average energy and fluctuations
     */
    void average(Energy &energy, Energy &fluctuation,
		 math::Matrix &pressure, 
		 math::Matrix &pressure_fluctuations);

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

    /**
     * the average pressure.
     */
    math::Matrix m_pressure_average;
    /**
     * the squared average pressure.
     */
    math::Matrix m_square_pressure_average;

  }; // Energy_Average
  
} // configuration

#endif

  
