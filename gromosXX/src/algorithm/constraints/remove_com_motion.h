/**
 * @file remove_com_motion.h
 * remove center of mass translational and angular momentum
 */

#ifndef INCLUDED_REMOVE_COM_MOTION_H
#define INCLUDED_REMOVE_COM_MOTION_H

namespace algorithm
{
  /**
   * @class Remove_COM_Motion
   * implements com removal
   */
  class Remove_COM_Motion : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Remove_COM_Motion(std::ostream & os = std::cout) : Algorithm("RemoveCOMMotion"), os(os) {}

    /**
     * Destructor.
     */
    virtual ~Remove_COM_Motion() {}
    
    /**
     * apply COM removal.
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
      os << "Remove centre of mass motion\n";
      return 0;
    };

  protected:
    std::ostream & os;
  };
  
} //algorithm

#endif
