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

      typedef interaction::Basic_Pairlist<t_simulation,
      interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, 
      interaction::Twinrange_Chargegroup_Filter<t_simulation,
      interaction::Nonbonded_Base,
      interaction::Nonbonded_Inner_Loop<
      t_simulation, interaction::Storage> > > > pairlist_type;
      
      typedef interaction::Nonbonded_Inner_Loop<t_simulation,
      typename t_simulation::system_type> innerloop_type;

      typedef interaction::Basic_Pairlist<t_simulation,
      interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, 
      interaction::Twinrange_Chargegroup_Filter<t_simulation,
      interaction::Nonbonded_Base,
      interaction::Nonbonded_Inner_Loop_Virial<
      t_simulation, interaction::Storage> > > > pairlist_virial_type;
      
      typedef interaction::Nonbonded_Inner_Loop_Virial<t_simulation,
      typename t_simulation::system_type> innerloop_virial_type;
    

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
     * the interactions:
     * quartic bond interaction.
     */
    interaction::Quartic_bond_interaction<t_simulation> * m_qbond_interaction;
    /**
     * angle interaction
     */
    interaction::angle_interaction<t_simulation> * m_angle_interaction;

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
