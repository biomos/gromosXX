/**
 * @file temperature/temperature_calculation.h
 * calculate the temperature.
 */

#ifndef INCLUDED_TEMPERATURE_CALCULATION_H
#define INCLUDED_TEMPERATURE_CALCULATION_H

namespace algorithm
{
  
  /**
   * @class Temperature_Calculation
   * temperature calculation.
   */
  class Temperature_Calculation : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Temperature_Calculation() : Algorithm("TemperatureCalculation") {}

    /**
     * Destructor.
     */
    virtual ~Temperature_Calculation() {}
    
    /**
     * apply the temperature calculation
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
		     bool quiet = false);
        
  private:

  };
  
} // algorithm

#endif
