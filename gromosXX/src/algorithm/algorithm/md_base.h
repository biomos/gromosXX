/**
 * @file md_base.h
 * the base for an md algorithm.
 */

#ifndef INCLUDED_MD_BASE_H
#define INCLUDED_MD_BASE_H

namespace algorithm
{
  /**
   * @class MD_Base
   * the base of an MD algorithm
   */
  template<typename t_md_spec=MD_spec<>,
	   typename t_interaction_spec=Interaction_spec<
    typename t_md_spec::simulation_type>
  >
  class MD_Base
  {
  public:
    /**
     * Constructor.
     */
    MD_Base();

    /**
     * Destructor.
     */
    virtual ~MD_Base();

    /**
     * simulation accessor.
     */
    typename t_md_spec::simulation_type & simulation();
    
    /**
     * forcefield accessor.
     */
    interaction::Forcefield<typename t_md_spec::simulation_type,
			    t_interaction_spec> & forcefield();
    
    /**
     * temperature coupling algorithm.
     */
    typename t_md_spec::temperature_type & temperature_algorithm();
    
    /**
     * pressure coupling algorithm.
     */
    typename t_md_spec::pressure_type & pressure_algorithm();
    
    /**
     * distance constraint algorithm.
     */
    typename t_md_spec::distance_constraint_type & distance_constraint_algorithm();

    /**
     * integration algorithm.
     */
    typename t_md_spec::integration_type & integration_algorithm();

    /**
     * the trajectory.
     */
    io::OutG96Trajectory<typename t_md_spec::simulation_type> & trajectory();

    /**
     * a title for the output files.
     */
    std::string title;

  protected:
    /**
     * trajecorty file
     */
    io::OutG96Trajectory<typename t_md_spec::simulation_type> *m_trajectory;
    /**
     * additional output file.
     */
    std::ostream * m_print_file;
    /**
     * the time step.
     */
    double m_dt;
    /**
     * simulation time.
     */
    double m_time;
    /**
     * print energy every .. steps.
     */
    int m_print_energy;
    /**
     * print pairlist every .. steps.
     */
    int m_print_pairlist;
    /**
     * print the force every .. steps.
     */
    int m_print_force;
    /**
     * remove com every .. steps.
     */
    int m_remove_com;
    /**
     * print com every .. steps.
     */
    int m_print_com;
    /**
     * calculate the pressure?
     * which kind of virial?
     */
    int m_calculate_pressure;      

    /**
     * are we performing a perturbation run.
     * which of course we cannot from this class
     * but (maybe???) from derived classes.
     */
    bool m_do_perturbation;

    /**
     * simulation.
     */
    typename t_md_spec::simulation_type m_simulation;
    /**
     * forcefield.
     */
    interaction::Forcefield<typename t_md_spec::simulation_type,
			    t_interaction_spec> m_forcefield;
    /**
     * temperature coupling.
     */
    typename t_md_spec::temperature_type m_temperature;
    /**
     * pressure coupling.
     */
    typename t_md_spec::pressure_type m_pressure;
    /**
     * distance constraint algorithm.
     */
    typename t_md_spec::distance_constraint_type m_distance_constraint;
      
    /**
     * integration algorithm.
     */
    typename t_md_spec::integration_type m_integration;

    /**
     * parse the print argument
     * for pairlist and force that
     * are not present in the input file...
     */
    virtual void parse_print_argument(io::Argument &args);

    /**
     * shake the initial positions and
     * velocities if required
     */
    void init_pos_vel(int init);
    
  };
  
} // algorithm

#include "md_base.tcc"

#endif

  
