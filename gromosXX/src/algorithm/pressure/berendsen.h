/**
 * @file pressure/berendsen.h
 * berendsen barostat
 */

#ifndef INCLUDED_PRESSURE_BERENDSEN_H
#define INCLUDED_PRESSURE_BERENDSEN_H

namespace algorithm
{
  
  /**
   * @class Berendsen_Barostat
   * the Berendsen barostat.
   */
  class Berendsen_Barostat
  {
  public:
    /**
     * Constructor.
     */
    Berendsen_Barostat();
    /**
     * apply the pressure scaling (weak coupling)
     * or just calculate the pressure (if tau=-1).
     * @param sim the simulation.
     * @param dt the time step.
     */
    template<typename t_simulation>
    void apply(t_simulation &sim, double const dt);

    /**
     * initialize.
     */
    void initialize(int ntp, math::Matrix pres0, double comp, double tau);
    
  private:

    int m_ntp;
    math::Matrix m_pres0;
    double m_comp;
    double m_tau;

  };
  
} // algorithm

#include "berendsen.tcc"

#endif
