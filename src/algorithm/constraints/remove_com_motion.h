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
   * implements centre of mass motion removal
   * input switches determine whether to remove
   * translational or angular centre of mass
   * momentum (or both).
   * For periodic boundary conditions, only
   * translational centre of mass motion is
   * removed.
   * It is recommended to remove it at every step.
   *
   * @TODO check if there is a Gromos96 bug in cenmas.
   * seems like absolute positions instead of
   * relative to centre of mass are taken in
   * angular momentum calculation.
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
		     bool quiet = false);

    /**
     * calculate and remove translational centre of mass motion
     */
    double remove_com_translation(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  bool remove_trans = true);
    
    /**
     * calculate and remove angular centre of mass motion
     */
    double remove_com_rotation(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       bool remove_rot = true);

    /**
     * add centre of mass rotation
     */
    double add_com_rotation(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    math::Vec com_L);
    

  protected:
    std::ostream & os;
  };
  
} //algorithm

#endif
