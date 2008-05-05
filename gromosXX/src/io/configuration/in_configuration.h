/**
 * @file in_configuration.h
 * read in a G96 trajectory file.
 */

#ifndef INCLUDED_IN_CONFIGURATION_H
#define INCLUDED_IN_CONFIGURATION_H

#include <gromosXX/io/configuration/inframe.h>

namespace util
{
  struct Replica_Data;
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
     * Read in replica exchange MD configurations
     */
    void read_replica
    (
     std::vector<configuration::Configuration> & conf, 
     topology::Topology &topo, 
     simulation::Simulation & sim,
     std::vector<util::Replica_Data> & replica_data,
     std::ostream & os = std::cout);

    /**
     * read the next frame, topology has to be already prepared.
     */
    bool read_next(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim,
		   std::ostream & os = std::cout);
    
    /**
     * read reference positions for position restraining
     * @param hasTitle the buffer has a title
     */
    static bool _read_refposition(
            std::vector<topology::position_restraint_struct> &res,
            std::vector<std::string> &buffer, bool hasTitle);
    /**
     * read B-factors for position restraining
     * @param hasTitle the buffer has a title
     */
    static bool _read_bfactor(
            std::vector<topology::position_restraint_struct> &res,
            std::vector<std::string> &buffer, bool hasTitle);

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
     * read distance restraint averages
     */
    bool read_distance_restraint_averages(topology::Topology &topo, 
                                          configuration::Configuration &conf, 
                                          simulation::Simulation & sim,
                                          std::ostream & os);
    
    /**
     * read position restraint averages and bfactors
     */
    bool read_position_restraints(topology::Topology &topo, 
                                  simulation::Simulation & sim,
                                  std::ostream & os);

    /**
     * read replica information
     */
    bool read_replica_information
    (
     std::vector<util::Replica_Data> & replica_data,
     std::ostream & os = std::cout
     );
				  
    /**
     * read POSITION block.
     */
    bool _read_position(math::VArray &pos, std::vector<std::string> &buffer,
			int const num);
    /**
     * read POSITIONRED block.
     */
    bool _read_positionred(math::VArray &pos, std::vector<std::string> &buffer,
			   int const num);
    
    /**
     * read COSPOSITION block.
     */
    bool _read_cos_position(math::VArray &pos, std::vector<std::string> &buffer,
			   int const num);
    
    /**
     * read VELOCITY block.
     */
    bool _read_velocity(math::VArray &vel, std::vector<std::string> &buffer,
			int const num);
    /**
     * read POSITIONRED block.
     */
    bool _read_velocityred(math::VArray &vel, std::vector<std::string> &buffer,
			   int const num);
    /**
     * read TRICLINICBOX block.
     */
    bool _read_box(math::Box &box, std::vector<std::string> &buffer,
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

    bool _read_stochastic_integral(math::VArray & stochastic_integral,
				   std::vector<std::string> &buffer,
				   int num_atoms, std::string &seed);
    
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
     * read time information.
     */
    bool _read_time(std::vector<std::string> &buffer,
		    double & t);
    /**
     * read time and step information.
     */
    bool _read_time_step(std::vector<std::string> &buffer,
			 simulation::Simulation & sim);
    
  };

} // io

#endif
