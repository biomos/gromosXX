/**
 * @file OutTrajectory.h
 * write a G96 trajectory file.
 */

#ifndef INCLUDED_OUTTRAJECTORY_H
#define INCLUDED_OUTTRAJECTORY_H

namespace io {
  
  /**
   * @enum format
   * output format.
   */
  enum output_format {
    /**
     * reduced format (the default).
     */
    reduced,
    /**
     * decorated format.
     */
    decorated,
    /**
     * final structure.
     */
    final
  };

  /**
   * @class OutTrajectory
   * writes out a trajectory file.
   */
  template<typename t_simulation>
  class OutTrajectory 
  {
  public:    
    /**
     * Constructor.
     */
    OutTrajectory(std::ostream &os, std::ostream &final);
    /**
     * write out the title block.
     */
    void print_title(std::string title);
    /**
     * write out a timestep.
     */
    OutTrajectory & operator<<(t_simulation &sim);
    /**
     * set output style.
     */
    OutTrajectory & operator<<(output_format f);
    /**
     * write a velocity trajectory.
     */
    void velocity_trajectory(std::ostream &os);
    /**
     * write a force trajectory.
     */
    void force_trajectory(std::ostream &os);
    
  private:
    output_format m_format;
    output_format m_old_format;
    std::ostream *m_pos_traj;
    std::ostream *m_final_traj;
    std::ostream *m_vel_traj;
    std::ostream *m_force_traj;
    
    bool m_pos;
    bool m_vel;
    bool m_force;

    void _print_timestep(t_simulation &sim, std::ostream &os);
    void _print_position(simulation::system &sys, simulation::topology &topo,
			 std::ostream &os);
    void _print_velocity(simulation::system &sys, simulation::topology &topo,
			 std::ostream &os);
    void _print_force(simulation::system &sys, simulation::topology &topo,
			 std::ostream &os);
    void _print_positionred(simulation::system &sys, std::ostream &os);
    void _print_velocityred(simulation::system &sys, std::ostream &os);
    void _print_forcered(simulation::system &sys, std::ostream &os);
    void _print_box(simulation::system &sys, std::ostream &os);
  };
  
} // io

// template methods
#include "OutTrajectory.tcc"

#endif

  
