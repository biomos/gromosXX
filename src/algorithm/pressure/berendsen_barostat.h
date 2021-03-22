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

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      // os << "Berendsen barostat\n";
      return 0;
    };

  private:

  };
  
} // algorithm

#endif
