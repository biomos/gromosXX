/**
 * @file OutTrajectory.h
 * write a G96 trajectory file.
 */

#ifndef INCLUDED_OUTTRAJECTORY_H
#define INCLUDED_OUTTRAJECTORY_H

namespace io {
  
  /**
   * @enum output_format
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
    OutTrajectory(std::ostream &os, std::ostream &final, 
		  int every=1, bool auto_delete = false);

    /**
     * Destructor.
     */
    ~OutTrajectory();
    
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
    void velocity_trajectory(std::ostream &os, int every=1);
    /**
     * write a force trajectory.
     */
    void force_trajectory(std::ostream &os, int every=1);
    /**
     * write an energy trajectory.
     */
    void energy_trajectory(std::ostream &os, int every=1);
    /**
     * precision of output.
     */
    void precision(int prec, int add=6);
    int precision();
    void force_precision(int prec, int add=9);
    int force_precision();
    
  protected:
    output_format m_format;
    output_format m_old_format;
    std::ostream *m_pos_traj;
    std::ostream *m_final_traj;
    std::ostream *m_vel_traj;
    std::ostream *m_force_traj;
    std::ostream *m_energy_traj;
    
    bool m_pos;
    bool m_vel;
    bool m_force;
    bool m_energy;

    int m_every_pos;
    int m_every_vel;
    int m_every_force;
    int m_every_energy;

    int m_precision;
    int m_force_precision;
    
    int m_width;
    int m_force_width;

    bool m_auto_delete;

    void _print_timestep(t_simulation &sim, std::ostream &os);

    void _print_old_timestep(t_simulation &sim, std::ostream &os);

    void _print_position(typename t_simulation::system_type &sys,
			 typename t_simulation::topology_type &topo,
			 std::ostream &os);

    void _print_velocity(typename t_simulation::system_type &sys,
			 typename t_simulation::topology_type &topo,
			 std::ostream &os);

    void _print_force(typename t_simulation::system_type &sys,
		      typename t_simulation::topology_type &topo,
		      std::ostream &os);

    void _print_positionred(typename t_simulation::system_type &sys, std::ostream &os);

    void _print_velocityred(typename t_simulation::system_type &sys, std::ostream &os);

    void _print_forcered(typename t_simulation::system_type &sys, std::ostream &os);
    void _print_energyred(t_simulation &sim, std::ostream &os);
    void _print_volumepressurered(typename t_simulation::system_type &sys, std::ostream &os);
    
    void _print_box(typename t_simulation::system_type &sys, std::ostream &os);
  };
  
} // io

// template methods
#include "OutTrajectory.tcc"

#endif

  
