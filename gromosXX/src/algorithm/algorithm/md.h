/**
 * @file md.h
 * the md algorithm.
 */

#ifndef INCLUDED_MD_H
#define INCLUDED_MD_H

namespace algorithm
{
  /**
   * @class MD
   * the MD algorithm
   */
  template<typename t_simulation,
	   typename t_temperature = algorithm::Berendsen_Thermostat,
	   typename t_pressure = algorithm::Berendsen_Barostat,
	   typename t_distance_constraint = algorithm::Shake<t_simulation>,
	   typename t_integration = algorithm::Leap_Frog<t_simulation> >
  class MD
  {
  public:
    typedef t_simulation simulation_type;
    typedef t_temperature temperature_algorithm_type;
    typedef t_pressure pressure_algortihm_type;
    typedef t_distance_constraint distance_constraint_type;
    typedef t_integration integration_algorithm_type;
    typedef interaction::twin_range_pairlist_cg<t_simulation> pairlist_type;
    
    /**
     * Constructor.
     */
    MD(t_simulation &sim);

    /**
     * Destructor.
     */
    virtual ~MD();

    /**
     * simulation accessor.
     */
    t_simulation & simulation();
    
    /**
     * forcefield accessor.
     */
    interaction::Forcefield<t_simulation> & forcefield();
    
    /**
     * temperature coupling algorithm.
     */
    t_temperature & temperature_algorithm();
    
    /**
     * pressure coupling algorithm.
     */
    t_pressure & pressure_algorithm();
    
    /**
     * distance constraint algorithm.
     */
    t_distance_constraint & distance_constraint_algorithm();

    /**
     * integration algorithm.
     */
    t_integration & integration_algorithm();

    /**
     * the trajectory.
     */
    io::OutG96Trajectory<t_simulation> & trajectory();
    
    /**
     * initialization.
     */
    int initialize(io::Argument &args);

    /**
     * initialize the perturbation parameters.
     */
    int init_perturbation(io::Argument &args);
    
    /**
     * run the system.
     * @param time the time to run the system.
     */
    void run(double time = -1);
    
    /**
     * create a Gromos96 like forcefield.
     */
    virtual void G96Forcefield(io::InTopology &topo,
			       io::InInput &input,
			       typename t_simulation::topology_type
			       &the_topology);


  protected:

    /**
     * simulation.
     */
    t_simulation & m_simulation;
    /**
     * forcefield.
     */
    interaction::Forcefield<t_simulation> m_forcefield;
    /**
     * temperature coupling.
     */
    t_temperature m_temperature;
    /**
     * pressure coupling.
     */
    t_pressure m_pressure;
    /**
     * distance constraint algorithm.
     */
    t_distance_constraint m_distance_constraint;
    /**
     * integration algorithm.
     */
    t_integration m_integration;
    
    /**
     * trajecorty file
     */
    io::OutG96Trajectory<t_simulation> *m_trajectory;
    /**
     * velocity trajectory
     */
    std::ofstream m_velocity_file;

    /**
     * force trajectory
     */
    std::ofstream m_force_file;
    
    /**
     * energy trajectory
     */
    std::ofstream m_energy_file;
  
    /**
     * trajectory file
     */
    std::ofstream m_trajectory_file;
    
    /**
     * final file
     */
    std::ofstream m_final_file;

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
     * calculate the pressure?
     * which kind of virial?
     */
    int m_calculate_pressure;
    
    /**
     * parse the print argument
     * for pairlist and force that
     * are not present in the input file...
     */
    void parse_print_argument(io::Argument &args);
    
  };
  
} // algorithm

#include "md.tcc"

#endif
