/**
 * @file in_configuration.h
 * read in a G96 trajectory file.
 */

#ifndef INCLUDED_IN_CONFIGURATION_H
#define INCLUDED_IN_CONFIGURATION_H

#include "inframe.h"

namespace util
{
  struct Umbrella;
  class BS_Umbrella;
}

namespace io {

  /**
   * @class In_Configuration
   * reads in a trajectory file and parses
   * it into configuration::Configuration
   */
  class In_Configuration : public In_Frame {  

  public:
    /**
     * Default Constructor.
     */
    In_Configuration() : In_Frame(){}

    /**
     * Constructor.
     */
    In_Configuration(std::istream& is) : In_Frame(is) {}

    /**
     * Read in a G96 trajectory into the Configuration.
     */
    void read(configuration::Configuration &conf, 
	      topology::Topology &topo, 
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

    /**
     * read the next frame, topology has to be already prepared.
     */
    bool read_next(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim,
		   std::ostream & os = std::cout,
                   bool do_read_time = true);
    
    /**
     * read B-factors for position restraining
     * @param hasTitle the buffer has a title
     */
    static bool _read_bfactor(
            math::SArray & b,
            std::vector<std::string> &buffer, bool hasTitle);

    /**
     * read local elevation umbrellas
     * @param reset reset the bias
     */
    static bool _read_leusbias(
            std::vector<util::Umbrella> & umbrellas,
            std::vector<std::string> &buffer, bool reset);

    /**
     * read POSITION block.
     */
    static bool _read_position(math::VArray &pos, std::vector<std::string> &buffer,
			topology::Topology &topo, configuration::Configuration & conf, std::string blockname = "POSITION");

    /**
     * Just try to read the positions without any checks.
     */
    bool read_position_plain(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       std::ostream & os = std::cout);

  private:
    /**
     * check coordinates / velocities
     */
    bool check_coordinates(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   int num_coords,
			   std::ostream & os = std::cout);

    /**
     * try to get positions
     */
    bool read_position(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       std::ostream & os = std::cout);
    /**
     * try to get cos positions
     */
    bool read_cos_position(topology::Topology & topo,
		           configuration::Configuration & conf,
		           simulation::Simulation & sim,
		           std::ostream & os = std::cout);
    
    /**
     * try to get velocities
     */
    bool read_velocity(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       std::ostream & os = std::cout);
    
    /**
     * try to get the lattice shifts
     */
    bool read_lattice_shifts(topology::Topology & topo,
		             configuration::Configuration & conf,
		             simulation::Simulation & sim,
		             std::ostream & os = std::cout);
    
    /**
     * try to get time, step
     */
    bool read_time(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim,
		   std::ostream & os = std::cout);

    /**
     * try to get time, step
     */
    bool read_time_step(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			std::ostream & os = std::cout);

    /**
     * try to get box
     */
    bool read_box(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  std::ostream & os = std::cout);

    /**
     * try to get jvalue average data
     */
    bool read_jvalue(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::ostream & os = std::cout);

    /**
     * try to get xray average data
     */
    bool read_xray(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::ostream & os = std::cout);

    /**
     * try to get pscale continuation data
     */
    bool read_pscale(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::ostream & os = std::cout);

    /**
     * try to get flexible constraints data
     */
    bool read_flexv(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::ostream & os = std::cout);
    
    /**
     * read stochastic integral
     */
    bool read_stochastic_integral(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  std::ostream & os = std::cout);
    /**
     * read perturbation restart data
     */
    bool read_perturbation(topology::Topology & topo,
            simulation::Simulation & sim,
            std::ostream & os = std::cout);
    /**
     * read distance restraint averages
     */
    bool read_distance_restraint_averages(topology::Topology &topo, 
                                          configuration::Configuration &conf, 
                                          simulation::Simulation & sim,
                                          std::ostream & os);
    /**
     * read Nose-Hoover-Chains variables
     */
    bool read_nose_hoover_chains(topology::Topology &topo, 
                                 configuration::Configuration &conf, 
                                 simulation::Simulation & sim,
                                 std::ostream & os);
    /**
     * read position restraint averages and bfactors
     */
    bool read_position_restraints(topology::Topology &topo,
                                  configuration::Configuration &conf, 
                                  simulation::Simulation & sim,
                                  std::ostream & os);

    /**
     * read configuration of roto translational constraints
     */
    bool read_rottrans(topology::Topology &topo, 
                       configuration::Configuration &conf, 
                       simulation::Simulation & sim,
                       std::ostream & os);

    /**
     * read order parameter restraint averages
     */
    bool read_order_parameter_restraint_averages(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation & sim,
            std::ostream & os);

    /**
     * read RDC average data
     */
    bool read_rdc(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::ostream & os = std::cout);

    /**
     * read LE umbrella potentials
     */
    bool read_leusbias(topology::Topology &topo,
                                  configuration::Configuration &conf,
                                  simulation::Simulation & sim,
                                  std::ostream & os);
    /**
     * Read the memory from the configuration file
     */
    bool read_bsleus(topology::Topology &topo,
                                  configuration::Configuration &conf,
                                  simulation::Simulation & sim,
                                  std::ostream & os);
		
    /**
     * read configuration of aeds parameter seach simulation
     */
    bool read_aedssearch(topology::Topology &topo,
                                  configuration::Configuration &conf,
                                  simulation::Simulation & sim,
                                  std::ostream & os);

     /**
      * read configuration of GAMD parameter search simulation
      */
     bool read_gamdstat(topology::Topology &topo,
                                  configuration::Configuration &conf,
                                  simulation::Simulation & sim,
                                  std::ostream & os);

    /**
     * Read in the memory function of the umbrellas (BSLEUSMEM block)
     * @param bs_umbrella
     * @param buffer
     * @return wheter successfull or not.
     */
    bool _read_bsleus(util::BS_Umbrella &bs_umbrella, 
                      std::vector<std::string> buffer);
    /**
     * Read in position in the BSLEUS subspace (BSLEUSPOS block)
     * @param bs_umbrella
     * @param buffer
     * @return wheter successfull or not.
     */
    bool _read_bsleuspos(util::BS_Umbrella &bs_umbrella, 
                      std::vector<std::string> buffer);
    /**
     * read POSITIONRED block.
     */
    bool _read_positionred(math::VArray &pos, std::vector<std::string> &buffer,
			   topology::Topology &topo, configuration::Configuration & conf);
    
    /**
     * read COSDISPLACEMENTS block.
     */
    bool _read_cos_position(math::VArray &pos, std::vector<std::string> &buffer,
			   int const num);
    
    /**
     * read VELOCITY block.
     */
    bool _read_velocity(math::VArray &vel, std::vector<std::string> &buffer,
			topology::Topology &topo);
    
    /**
     * read LATTICESHIFTS block.
     */
    bool _read_lattice_shifts(math::VArray &vel, std::vector<std::string> &buffer,
			topology::Topology &topo);
    
    /**
     * read POSITIONRED block.
     */
    bool _read_velocityred(math::VArray &vel, std::vector<std::string> &buffer,
			   topology::Topology &topo);
    /**
     * read TRICLINICBOX block.
     */
    bool _read_box(math::Box &box, double &phi, double &theta, double &psi, 
                   std::vector<std::string> &buffer,
		   math::boundary_enum const boundary);
    /**
     * read GENBOX block.
     */
    bool _read_genbox(math::Box &box, double &phi, double &theta, double &psi, 
                   std::vector<std::string> &buffer,
		   math::boundary_enum const boundary);
    /**
     * read BOX block.
     */
    bool _read_g96_box(math::Box &box, std::vector<std::string> &buffer);

    /**
     * read FLEXV block.
     */
    bool _read_flexv(std::vector<double> &flexv,
		     std::vector<std::string> &buffer,
		     std::vector<topology::two_body_term_struct>
		     const & constr,
		     std::vector<topology::perturbed_two_body_term_struct>
		     const & pert_constr);

    /**
     * read STOCHINT block
     */
    bool _read_stochastic_integral(math::VArray & stochastic_integral,
				   std::vector<std::string> &buffer,
				   int num_atoms, std::string &seed);
    /**
     * read PERTDATA block
     */
    bool _read_pertdata(topology::Topology & topo,
                        std::vector<std::string> & buffer);
    
    /**
     * read DISRESEXPAVE block
     */
    bool _read_distance_restraint_averages(
                 std::vector<std::string> &buffer,
                 const std::vector<topology::distance_restraint_struct> & distanceress,
                 std::vector<double> &distanceres_av);
    
    /**
     * read JVALUE averages.
     */
    bool _read_jvalue_av(std::vector<std::string> &buffer,
			 std::vector<double> & jvalue_av,
			 std::vector<topology::jvalue_restraint_struct> const & jval_res);

    /**
     * read XRAY averages.
     */
    bool _read_xray_av(std::vector<std::string> &buffer,
			 std::vector<configuration::Configuration::special_struct::xray_struct> & xray_av,
			 std::vector<topology::xray_restraint_struct> const & xray_res,
                         std::vector<topology::xray_restraint_struct> const & xray_rfree);

    /**
     * read X-ray umbrella weight threshold data
     */
    bool _read_xray_umbrellaweightthesholds(std::vector<std::string> &buffer,
                                            std::vector<topology::xray_umbrella_weight_struct> & umb_weight);

    /**
     * read the B factors/occupancies from the configuration
     */
    bool _read_xray_bfactors(std::vector<std::string> &buffer,
            std::vector<configuration::Configuration::special_struct::xray_bfoc_struct> & bfoc);

    /**
     * read JVALUE local elevation epsilons.
     */
    bool _read_jvalue_le(std::vector<std::string> &buffer,
			 std::vector<std::vector<double> > & jvalue_epsilon,
			 std::vector<topology::jvalue_restraint_struct> const & jval_res,
                         unsigned int const & grid_size);
    
    /**
     * read Periodic Scaling (PSCALE) data.
     */
    bool _read_pscale_jrest(std::vector<std::string> &buffer,
			    configuration::Configuration::special_struct::pscale_struct &pscale,
			    std::vector<topology::jvalue_restraint_struct> const & jval_res);

    /**
     * read ORDERPARAMRESEXPAVE block
     */
    bool _read_order_parameter_restraint_averages(
                 std::vector<std::string> & buffer,
                 const std::vector<topology::order_parameter_restraint_struct> & oparamres,
                 std::vector<math::Matrix> & Q_avg,
                 std::vector<double> & D_avg);
    /**
     * read ORDERPARAMRESWINAVE block
     */
    bool _read_order_parameter_restraint_average_window(
                 std::vector<std::string> & buffer,
                 unsigned int window_size,
                 const std::vector<topology::order_parameter_restraint_struct> & oparamres,
                 std::vector<std::list<math::Matrix> > & Q_avg,
                 std::vector<std::list<double> > & D_avg);


    /**
     * read RDC averages
     */
    bool _read_rdc_av(std::vector<std::string> &buffer,
             std::vector<configuration::Configuration::special_struct::rdc_struct> &rdc,
			 std::vector<std::vector<topology::rdc_restraint_struct> > const &rdc_res,
		     std::ostream & os = std::cout);

    /**
     * read RDC magnetic field vectors
     */
    bool _read_rdc_mf(std::vector<std::string> &buffer,
             std::vector<configuration::Configuration::special_struct::rdc_struct> &rdc,
			 std::vector<std::vector<topology::rdc_restraint_struct> > const &rdc_res,
		     std::ostream & os = std::cout);
    
    /**
     * read RDC alignment tensor
     */
    bool _read_rdc_t(std::vector<std::string> &buffer,
             std::vector<configuration::Configuration::special_struct::rdc_struct> &rdc,
		     std::ostream & os = std::cout);
    
    /**
     * read RDC spherical harmonics
     */
    bool _read_rdc_sh(std::vector<std::string> &buffer,
             std::vector<configuration::Configuration::special_struct::rdc_struct> &rdc,
		     std::ostream & os = std::cout);

    /**
     * read RDC stochastic integrals
     */
    bool _read_rdc_stochint(std::vector<std::string> &buffer,
             std::vector<configuration::Configuration::special_struct::rdc_struct> &rdc,
             simulation::rdc_type_enum & mode,
		     std::ostream & os = std::cout);

    /**
     * read time information
     */
    bool _read_time(std::vector<std::string> &buffer,
		    double & t);
    /**
     * read time and step information.
     */
    bool _read_time_step(std::vector<std::string> &buffer,
			 simulation::Simulation & sim);
    /**
     * read NHCVARIABLES block
     */
    bool _read_nose_hoover_chain_variables(std::vector<std::string> &buffer,
            simulation::Multibath & multibath);
    
    /**
     * read ROTTRANSREFPOS block
     */
    bool _read_rottrans(std::vector<std::string> &buffer, unsigned int last,
            configuration::Configuration::special_struct::rottrans_constr_struct & rottrans);
    
    /**
     * read AEDSSEARCH block
     */
    bool _read_aedssearch(std::vector<std::string> &buffer,
            simulation::Simulation & sim, unsigned int numstates);

    /**
     * read GAMDSTAT block
     */
    bool _read_gamdstat(std::vector<std::string> &buffer,
            simulation::Simulation & sim, unsigned int numagroups);
  };

} // io

#endif
