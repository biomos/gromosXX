/**
 * @file in_configuration.h
 * read in a G96 trajectory file.
 */

#ifndef INCLUDED_IN_CONFIGURATION_H
#define INCLUDED_IN_CONFIGURATION_H

#include <gromosXX/io/configuration/inframe.h>

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
    void read(configuration::Configuration &conf, topology::Topology &topo, 
	      simulation::Parameter const &param);
  private:
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
    
    /**
     * read JVALUE averages.
     */
    bool _read_jvalue_av(std::vector<std::string> &buffer,
			 std::vector<double> & jvalue_av,
			 std::vector<topology::jvalue_restraint_struct> const & jval_res);

    /**
     * read Periodic Scaling (PSCALE) data.
     */
    bool _read_pscale_jrest(std::vector<std::string> &buffer,
			    configuration::Configuration::special_struct::pscale_struct &pscale,
			    std::vector<topology::jvalue_restraint_struct> const & jval_res);

  };
  

} // io

#endif
