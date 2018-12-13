/**
 * @file lattice_shift.h
 * keeping track of lattice shifts
 */

#ifndef INCLUDED_LATTICE_SHIFT_H
#define INCLUDED_LATTICE_SHIFT_H

namespace algorithm
{
  /**
   * @class Lattice_Shift_Tracker
   * keeps track of lattice shifts
   */
  class Lattice_Shift_Tracker : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Lattice_Shift_Tracker() : Algorithm("Lattice_Shift_Tracker") {}

    /**
     * Destructor.
     */
    virtual ~Lattice_Shift_Tracker(){}
    
    /**
     * put CG into box and keep track of shift
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
    
  protected:
    template<math::boundary_enum b>
    void _apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);
  };
}
#endif	/* INCLUDED_LATTICE_SHIFT_H */

