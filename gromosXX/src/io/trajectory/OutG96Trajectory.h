/**
 * @file OutG96Trajectory.h
 * write out a trajectory file,
 * gathering the positions according to
 * gromos 96 with respect to chargegroup center
 * of geometries...
 * (not all atoms independently)
 */

#ifndef INCLUDED_OUTG96TRAJECTORY_H
#define INCLUDED_OUTG96TRAJECTORY_H

namespace io
{
  
  /**
   * @class OutG96Trajectory
   * writes a trajectory in gromos 96 style.
   * The difference to standard writing of a
   * trajectory is, that we need to put the
   * chargegroup center of geometries into
   * the central box, not the atoms independently.
   */
  template<typename t_simulation>
  class OutG96Trajectory : public OutTrajectory<t_simulation>
  {
  public:
    /**
     * Constructor.
     */
    OutG96Trajectory(std::ostream &os, std::ostream &final, int every=1, 
		     bool auto_delete = false);
    /**
     * write out a timestep.
     */
    OutG96Trajectory & operator<<(t_simulation &sim);
    /**
     * set output style.
     */
    OutG96Trajectory & operator<<(output_format f);

    
  protected:

    void _print_position(typename t_simulation::system_type &sys,
			 typename t_simulation::topology_type &topo,
			 std::ostream &os);

    void _print_positionred(typename t_simulation::system_type &sys,
			    typename t_simulation::topology_type &topo,
			    std::ostream &os);
  };
  
} // io

// template methods
#include "OutG96Trajectory.tcc"

#endif
  
