/**
 * @file out_configuration.h
 * write a G96 trajectory file.
 */
/**
 * @page output output files
 * @date 28-10-2008
 * @section energy_trajectory energy trajectory
 *  - @ref energyredhelper
 *  - @ref volumepressurered
 * @section free_energy_trajectory free energy trajectory
 *  - see @ref energy_trajectory
 * @section block_averaged_energy_trajectory block averaged energy trajectory
 *  - see @ref energy_trajectory
 * @section free_energy_trajectory free energy trajectory
 *  - see @ref energy_trajectory
 * @section block_averaged_free_energy_trajectory block averaged free energy trajectory
 *  - see @ref energy_trajectory
 */
#ifndef INCLUDED_OUT_CONFIGURATION_H
#define INCLUDED_OUT_CONFIGURATION_H

namespace util {
  struct Replica_Data;
}

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
  class Out_Configuration {
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
    void trajectory(std::string const name, int every = 1);
    /**
     * write a velocity trajectory.
     */
    void velocity_trajectory(std::string const name, int every = 1);
    /**
     * write a force trajectory.
     */
    void force_trajectory(std::string const name, int every = 1);
    /**
     * write a special trajectory.
     */
    void special_trajectory(std::string const name, int every_cos = 1,
            int every_jvalue = 1, int every_xray = 1, int every_disres = 1, int every_disfieldres = 1,
            int every_angres = 1,
            int every_dihres = 1, int every_dat = 1, int every_leus = 1, int every_dipole = 1, 
            int every_current = 1, int every_adde = 1, int every_nemd = 1,
	    int every_oparam = 1, int every_rdc = 1, int every_bsleus = 1);
    /**
     * write an energy trajectory.
     */
    void energy_trajectory(std::string const name, int every = 1);
    /**
     * write a free energy trajectory
     */
    void free_energy_trajectory(std::string const name, int every = 1);
    /**
     * write block averaged energy / pressure / volume properties.
     */
    void block_averaged_energy(std::string const name, int every = 1);
    /**
     * write block averaged free energies.
     */
    void block_averaged_free_energy(std::string const name, int every = 1);
    /**
     * precision of output.
     */
    void precision(int prec, int add = 6);
    /**
     * precision accessor.
     */
    int precision();
    /**
     * force write precision.
     */
    void force_precision(int prec, int add = 9);
    /**
     * force write precision accessor.
     */
    int force_precision();

    /**
     * get the title.
     */
    std::string title() const {
      return m_title;
    }

    /**
     * set the title.
     */
    void title(std::string s) {
      m_title = s;
    }

    /**
     * get standard output stream.
     */
    std::ostream & output() {
      return m_output;
    }

    /**
     * print timestep.
     */
    void print_timestep(simulation::Simulation const &sim,
            std::ostream &os) {
      _print_timestep(sim, os);
    }
    
    /**
     * print timestep to special traj if any special terms are to be written
     */
    void _print_special_timestep(simulation::Simulation const &sim);

    // make them available for scripting!
    void _print_title(std::string title, std::string name,
            std::ostream &os);
    
    /*
     * Prints the ENEVERSION block for the (free) energy trajectories
     */
    void _print_ene_version(std::ostream &os);

    void _print_timestep(simulation::Simulation const &sim,
            std::ostream &os);

    void _print_old_timestep(simulation::Simulation const &sim,
            std::ostream &os);

    void _print_position(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_shake_failure(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_velocity(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_lattice_shifts(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_force(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os,
            bool constraint_force);

    void _print_positionred(configuration::Configuration const &conf,
            topology::Topology const &topo,
            int num,
            std::ostream &os);

    void _print_cos_position(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_velocityred(configuration::Configuration const &conf,
            topology::Topology const &topo,
            int num,
            std::ostream &os);

    void _print_forcered(configuration::Configuration const &conf,
            topology::Topology const &topo,
            int num,
            std::ostream &os,
            bool constraint_force);

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

    void _print_stochastic_integral(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_dihangle_trans(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_aedssearch(configuration::Configuration const &conf,
            simulation::Simulation const &sim,
            std::ostream &os);

     // ORIOL_GAMD
     void _print_gamdstat(configuration::Configuration const &conf,
            simulation::Simulation const &sim,
            std::ostream &os);

    /**
     * perturbation (slow groth) restart data
     */
    void _print_pertdata(topology::Topology const &topo,
            std::ostream &os);

    void _print_distance_restraints(
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_distance_restraint_averages(
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_disfield_restraints(
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_disfield_grid(
	    simulation::Parameter const &param,
	    configuration::Configuration const &conf,
	    topology::Topology const &topo,
	    std::ostream &os);

    void _print_angle_restraints(
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_dihedral_restraints(
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);
    
    void _print_position_restraints(simulation::Simulation const &sim,
            topology::Topology const &topo,
            configuration::Configuration const &conf,
            std::ostream &os);

    void _print_jvalue(simulation::Parameter const & param,
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os,
            bool formated = false);

    void _print_xray(simulation::Parameter const & param,
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os,
            bool final = false);

    void _print_xray_rvalue(simulation::Parameter const & param,
            configuration::Configuration const &conf,
            std::ostream &os);

    void _print_xray_umbrellaweightthresholds(simulation::Parameter const & param,
        topology::Topology const & topo,
        std::ostream &os);

    void _print_xray_bfactors(simulation::Parameter const & param,
        configuration::Configuration const &conf,
        std::ostream &os);

    void _print_pscale_jrest(configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_order_parameter_restraints(
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_order_parameter_restraint_averages(
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_order_parameter_restraint_average_window(
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os);

    void _print_rdc_values(simulation::Parameter const & param,
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os,
            bool formatted = false);

    void _print_rdc_averages(simulation::Parameter const & param,
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os,
            bool formatted = false);

    void _print_rdc_representation(simulation::Parameter const & param,
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os,
            bool formatted = false);

    void _print_rdc_stochastic_integrals(simulation::Parameter const & param,
            configuration::Configuration const &conf,
            topology::Topology const &topo,
            std::ostream &os,
            bool formatted = false);

    void _print_blockaveraged_energyred(configuration::Configuration const & conf,
            std::ostream & os);

    void _print_blockaveraged_volumepressurered(configuration::Configuration const & conf,
            simulation::Simulation const & sim,
            std::ostream & os);

    void _print_blockaveraged_free_energyred(configuration::Configuration const & conf,
            double dlamt,
            std::ostream & os);

    void _print_nose_hoover_chain_variables(const simulation::Multibath & multibath,
            std::ostream & os);

    void _print_rottrans(configuration::Configuration const &conf,
            simulation::Simulation const &sim,
            std::ostream &os);

    void _print_umbrellas(configuration::Configuration const & conf,
            std::ostream & os);
    
    void _print_bsleusmem(configuration::Configuration const &conf,
            std::ostream &os);
    
    void _print_bsleuspos(configuration::Configuration const &conf,
            std::ostream &os);
    
    void _print_bsleus(configuration::Configuration const &conf,
            std::ostream &os);

    template<math::boundary_enum b>
    void _print_dipole(simulation::Simulation const & sim,
                       topology::Topology const &topo,
                       configuration::Configuration const & conf,
                       std::ostream & os);

    void _print_current(simulation::Simulation const & sim,
                        topology::Topology const &topo,
                        configuration::Configuration const & conf,
                        std::ostream & os);
    
    void _print_adde(simulation::Simulation const & sim,
                        topology::Topology const &topo,
                        configuration::Configuration const & conf,
                        std::ostream & os);
    
    void _print_nemd(simulation::Simulation const & sim,
                        topology::Topology const &topo,
                        configuration::Configuration const & conf,
                        std::ostream & os);


  protected:
    std::ofstream m_pos_traj;
    std::ofstream m_final_conf;
    std::ofstream m_vel_traj;
    std::ofstream m_force_traj;
    std::ofstream m_energy_traj;
    bool m_has_replica_traj;
    std::ofstream m_replica_traj;
    std::ofstream m_free_energy_traj;
    std::ofstream m_blockaveraged_energy;
    std::ofstream m_blockaveraged_free_energy;
    std::ofstream m_special_traj;

    std::ostream & m_output;

    bool m_final;
    bool m_replica;

    /**
     * true if a special trajectory has to be written
     */
    bool m_write_special;

    int m_every_pos;
    int m_every_vel;
    int m_every_force;
    int m_every_energy;
    int m_every_free_energy;
    int m_every_blockaverage;
    
    int m_every_cos_pos;
    int m_every_jvalue;
    int m_every_xray;
    int m_every_disres;
    int m_every_disfieldres;
    int m_every_angres;
    int m_every_dihres;
    int m_every_dat;
    int m_every_leus;
    int m_every_bsleus;
    int m_every_dipole;
    int m_every_current;
    int m_every_adde;
    int m_every_nemd;
    int m_every_oparam;
    int m_every_rdc;

    bool m_write_blockaverage_energy;
    bool m_write_blockaverage_free_energy;

    int m_precision;
    int m_force_precision;
    int m_distance_restraint_precision;
    int m_disfield_restraint_precision;
    int m_angle_restraint_precision;
    int m_dihedral_restraint_precision;

    int m_width;
    int m_force_width;

    std::string m_title;
    /**
     * minimum energy for NTWSE trajectory
     */
    double minimum_energy;
    /**
     * The version of the energy trajectory, will be printed in the ENEVERSION
     * block. ene_ana compares the ENEVERSION of the trajectory with the one of
     * the ene_ana library used.
     * If modifications in the energy trajectory are done, this version number
     * should be adapted, along with an updated version of the ene_ana library.
     * The standard for the version number is the current date, in the format
     * YYYY-MM-DD, followed by a dash and an (incrementation) letter, for the 
     * case that two or more version should be comitted ona specific day.
     * Example:
     *   2015-06-23-A
     * (As discussed at GROMOS Meeting 22. June 2015)
     * The string is initialized at the head of the out_configuration.cc file.
     */
    static const std::string ene_version;
  };

} // io

#endif


