/**
 * @file energy_calculation.h
 * calculates the (total) energies and updates the averages
 */

#ifndef INCLUDED_ENERGY_CALCULATION_H
#define INCLUDED_ENERGY_CALCULATION_H

namespace algorithm
{
  /**
   * @class Energy_Calculation
   * calculates total energies, updates the averages
   */
  class Energy_Calculation : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Energy_Calculation() : Algorithm("EnergyCalculation") {}

    /**
     * Destructor.
     */
    virtual ~Energy_Calculation(){}
    
    /**
     * calculate the totals
     * update the averages
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    

  };
  
} // algorithm

#endif

