/**
 * @file pressure/berendsen_barostat.h
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
  class Berendsen_Barostat : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Berendsen_Barostat() : Algorithm("BerendsenBarostat") {}
    /**
     * Destructor.
     */
    virtual ~Berendsen_Barostat() {}
    
    /**
     * apply the pressure scaling (weak coupling)
     * or just calculate the pressure (if tau=-1).
     * @param sim the simulation.
     * @param dt the time step.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

  private:

  };
  
} // algorithm

#endif
