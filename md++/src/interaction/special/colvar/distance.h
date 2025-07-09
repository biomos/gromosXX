/**
 * @file distance.h
 * contact number restraining
 */

#ifndef INCLUDED_DISTANCE_RESTRAINT_INTERACTION_H
#define INCLUDED_DISTANCE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class Distance_Colvar
   * calculates distance and derivatives with respect to the position
   * which can then be used for energy and force calculation in the colvar
   * restraint interaction
   */
  class Distance_Colvar : public Colvar
  {
  public:
    /**
     * Constructor.
     */
    Distance_Colvar() : Colvar("Distance") {}
    
    /**
     * Destructor.
     */
    virtual ~Distance_Colvar() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
    /**
     * calculate distance and derivatives
     */
    virtual int calculate(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
				       
    topology::distance_restraint_struct_colvar *params;
    
  };
  
} // interaction

#endif
