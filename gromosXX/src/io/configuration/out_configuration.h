/**
 * @file out_configuration.h
 * write a G96 trajectory file.
 */

#ifndef INCLUDED_OUT_CONFIGURATION_H
#define INCLUDED_OUT_CONFIGURATION_H

namespace io {
  
  class Argument;

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
   * @class Out_Configuration
   * writes out the configuration to a file.
   */
  class Out_Configuration
  {
  public:    
    /**
     * Constructor.
     */
    Out_Configuration(std::string title,
		      std::ostream & os = std::cout);

    /**
     * Destructor.
     */
    ~Out_Configuration();
    
    /**
     * initialise the files.
     */
    void init(io::Argument & args, simulation::Parameter const & param);
    
    /**
     * write out a timestep.
     */
    void write(configuration::Configuration & conf,
	       topology::Topology const & topo,
	       simulation::Simulation const &sim,
	       output_format const form = reduced);

    /**
     * print out data (per time step).
     */
    void print(topology::Topology const & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation const & sim);

    /**
     * print out final data
     */
    void print_final(topology::Topology const & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation const & sim);

    /**
     * final structure.
     */
    void final_configuration(std::string const name);
    /**
     * write a trajectory.
     */
    void trajectory(std::string const name, int every=1);
    /**
     * write a velocity trajectory.
     */
    void velocity_trajectory(std::string const name, int every=1);
    /**
     * write a force trajectory.
     */
    void force_trajectory(std::string const name, int every=1);
    /**
     * write an energy trajectory.
     */
    void energy_trajectory(std::string const name, int every=1);
    /**
     * write a free energy trajectory
     */
    void free_energy_trajectory(std::string const name, int every=1);
    /**
     * write block averaged energy / pressure / volume properties.
     */
    void block_averaged_energy(std::string const name, int every=1);
    /**
     * write block averaged free energies.
     */
    void block_averaged_free_energy(std::string const name, int every=1);
    /**
     * precision of output.
     */
    void precision(int prec, int add=6);
    /**
     * precision accessor.
     */
    int precision();
    /**
     * force write precision.
     */
    void force_precision(int prec, int add=9);
    /**
     * force write precision accessor.
     */
    int force_precision();
    /**
     * get the title.
     */
    std::string title() const { return m_title; }
    /**
     * set the title.
     */
    void title(std::string s) { m_title = s; }
    /**
     * get standard output stream.
     */
    std::ostream & output() { return m_output; }
    /**
     * print timestep.
     */
    void print_timestep(simulation::Simulation const &sim,
			std::ostream &os)
    {
      _print_timestep(sim, os);
    }

    // make them available for scripting!
    void _print_title(std::string title, std::string name, 
		      std::ostream &os);

    void _print_timestep(simulation::Simulation const &sim,
			 std::ostream &os);

    void _print_old_timestep(simulation::Simulation const &sim,
			     std::ostream &os);

    void _print_position(configuration::Configuration const &conf,
			 topology::Topology const &topo,
			 std::ostream &os);

    void _print_velocity(configuration::Configuration const &conf,
			 topology::Topology const &topo,
			 std::ostream &os);

    void _print_force(configuration::Configuration const &conf,
		      topology::Topology const &topo,
		      std::ostream &os);

    void _print_positionred(configuration::Configuration const &conf,
			    topology::Topology const &topo,
			    std::ostream &os);

    void _print_velocityred(configuration::Configuration const &conf, 
			    std::ostream &os);

    void _print_forcered(configuration::Configuration const &conf, 
			 std::ostream &os);

    void _print_energyred(configuration::Configuration const &conf,
			  std::ostream &os);

    void _print_volumepressurered(topology::Topology const &topo,
				  configuration::Configuration const &conf, 
				  simulation::Simulation const &sim,
				  std::ostream &os);

    void _print_free_energyred(configuration::Configuration const &conf,
			       topology::Topology const &topo,
			       std::ostream &os);
    
    void _print_box(configuration::Configuration const &conf, 
		    std::ostream &os);

    void _print_flexv(configuration::Configuration const &conf,
		      topology::Topology const &topo,
		      std::ostream &os);

    void _print_jvalue(configuration::Configuration const &conf,
		       topology::Topology const &topo,
		       std::ostream &os);

    void _print_pscale_jrest(configuration::Configuration const &conf,
			     topology::Topology const &topo,
			     std::ostream &os);

    void _print_blockaveraged_energyred(configuration::Configuration const & conf,
					std::ostream & os);

    void _print_blockaveraged_volumepressurered(configuration::Configuration const & conf, 
						simulation::Simulation const & sim,
						std::ostream & os);

    void _print_blockaveraged_free_energyred(configuration::Configuration const & conf,
					     double dlamt,
					     std::ostream & os);

    void _print_replica_data(configuration::Configuration const & conf,
			     simulation::Simulation const & sim,
			     std::ostream & os);

  protected:
    std::ofstream m_pos_traj;
    std::ofstream m_final_conf;
    std::ofstream m_vel_traj;
    std::ofstream m_force_traj;
    std::ofstream m_energy_traj;
    std::ofstream m_free_energy_traj;
    std::ofstream m_blockaveraged_energy;
    std::ofstream m_blockaveraged_free_energy;
    std::ofstream m_replica_data;
    
    std::ostream & m_output;
    
    bool m_final;
    bool m_replica;
    
    int m_every_pos;
    int m_every_vel;
    int m_every_force;
    int m_every_energy;
    int m_every_free_energy;
    int m_every_blockaverage;

    bool m_write_blockaverage_energy;
    bool m_write_blockaverage_free_energy;

    int m_precision;
    int m_force_precision;
    
    int m_width;
    int m_force_width;

    std::string m_title;


  };
  
} // io

#endif

  
