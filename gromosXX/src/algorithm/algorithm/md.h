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
  template<typename t_spec=MD_spec>
  class MD
  {
  public:
    /**
     * Constructor.
     */
    MD();

    /**
     * Destructor.
     */
    virtual ~MD();

    /**
     * simulation accessor.
     */
    typename t_spec::simulation_type & simulation();
    
    /**
     * forcefield accessor.
     */
    interaction::Forcefield<typename t_spec::simulation_type> & forcefield();
    
    /**
     * temperature coupling algorithm.
     */
    typename t_spec::temperature_type & temperature_algorithm();
    
    /**
     * pressure coupling algorithm.
     */
    typename t_spec::pressure_type & pressure_algorithm();
    
    /**
     * distance constraint algorithm.
     */
    typename t_spec::distance_constraint_type & distance_constraint_algorithm();

    /**
     * integration algorithm.
     */
    typename t_spec::integration_type & integration_algorithm();

    /**
     * the trajectory.
     */
    io::OutG96Trajectory<typename t_spec::simulation_type> & trajectory();
    
    /**
     * initialization.
     */
    int initialize(io::Argument &args);
    
    /**
     * perform an md simulation.
     * calls run.
     */
    int do_md(io::Argument &args);

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
			       io::Argument &args);


    std::string title;
      
  protected:

    /**
     * simulation.
     */
    typename t_spec::simulation_type m_simulation;
    /**
     * forcefield.
     */
    interaction::Forcefield<typename t_spec::simulation_type> m_forcefield;
    /**
     * temperature coupling.
     */
    typename t_spec::temperature_type m_temperature;
    /**
     * pressure coupling.
     */
    typename t_spec::pressure_type m_pressure;
    /**
     * distance constraint algorithm.
     */
    typename t_spec::distance_constraint_type m_distance_constraint;
      
    /**
     * integration algorithm.
     */
    typename t_spec::integration_type m_integration;
    /**
     * trajecorty file
     */
    io::OutG96Trajectory<typename t_spec::simulation_type> *m_trajectory;
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
     * print com every .. steps.
     */
    int m_print_com;
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
     * are we performing a perturbation run.
     * which of course we cannot from this class
     * but (maybe???) from derived classes.
     */
    bool m_do_perturbation;
    
    /**
     * parse the print argument
     * for pairlist and force that
     * are not present in the input file...
     */
    void parse_print_argument(io::Argument &args);
    /**
     * open the input files.
     */
    void open_files(io::Argument &args, io::InTopology &topo,
		    io::InTrajectory &sys, io::InInput &input);
    /**
     * initialize the input.
     */
    virtual void init_input(io::Argument &args, io::InTopology &topo,
			    io::InTrajectory &sys, io::InInput &input);
    /**
     * read the input and setup a standard simulation.
     */
    virtual void read_input(io::Argument &args, io::InTopology &topo,
			    io::InTrajectory &sys, io::InInput &input);
    /**
     * initialize the output.
     */
    virtual void init_output(io::Argument &args, io::InInput &input);
    /**
     * print pairlists.
     */
    virtual void print_pairlists();
    /**
     * calculate and print the energies.
     */
    virtual void do_energies();
    
  };
  
} // algorithm

#include "md.tcc"

#endif
