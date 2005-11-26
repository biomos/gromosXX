/**
 * @file pressure_calculation.h
 * pressure calculation
 */

#ifndef INCLUDED_PRESSURE_CALCULATION_H
#define INCLUDED_PRESSURE_CALCULATION_H

namespace algorithm
{
  
  /**
   * @class Pressure_Calculation
   * pressure calculation.
   */
  class Pressure_Calculation : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Pressure_Calculation() : Algorithm("PressureCalculation") {}

    /**
     * Destructor.
     */
    virtual ~Pressure_Calculation() {}
    
    /**
     * apply the pressure calculation
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
      // os << "Pressure calculation\n";
      return 0;
    };

  private:

  };
  
} // algorithm

#endif
