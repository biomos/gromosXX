/**
 * @file berendsen_barostat.h
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
     * apply weak coupling.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

  private:

  };
  
} // algorithm

#endif
