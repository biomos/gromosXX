/**
 * @file print_block.h
 * routines to print out the various blocks.
 */

#ifndef INCLUDED_PRINT_BLOCK_H
#define INCLUDED_PRINT_BLOCK_H

namespace util
{
  struct Replica_Data;
}

namespace io
{
  /**
   * Print the dof coupling table from the MULTIBATH block.
   */
  void print_MULTIBATH_COUPLING(std::ostream &os, simulation::Multibath const &bath);

  /** 
   * Print the MULTIBATH block.
   */
  void print_MULTIBATH(std::ostream &os, simulation::Multibath const &bath,
		       configuration::Energy const &energy,
		       std::string title="TEMPERATURES");

  /** 
   * Print the DEGREES OF FREEDOM block.
   */
  void print_DEGREESOFFREEDOM(std::ostream &os, simulation::Multibath const &bath);

  /**
   * Print the PCOUPLE block.
   */
  void print_PCOUPLE(std::ostream &os, bool calc, 
		     math::pressure_scale_enum scale,
		     math::Matrix pres0, double comp, 
		     double tau, math::virial_enum vir,
                     int x_semi, int y_semi, int z_semi);

  /**
   * Print the PRESSURE block.
   */
  void print_PRESSURE(std::ostream &os, configuration::Configuration const & conf);
  
  /**
   * Print the ENERGY block.
   */
  void print_ENERGY(std::ostream &os, configuration::Energy const &e,
		    std::vector<unsigned int> const & energy_groups,
		    std::string const title = "ENERGIES",
		    std::string const type = "E_");
  
  /**
   * Print a matrix.
   */
  void print_MATRIX(std::ostream &os, math::Matrix const &m,
		    std::string const title);

  /**
   * Print the TIMESTEP block.
   */
  void print_TIMESTEP(std::ostream &os, double const steps, double const time);

  /**
   * Print the CENTREOFMASS block.
   */
  void print_CENTREOFMASS(std::ostream &os, double const ekin_trans, double const ekin_rot);
  
  /**
   * print replica exchange information
   * will precede standard blocks in trajectories,
   * identifying the current replica
   * reeval 0: data is from Ti,li
   *        1: data is from Tj,lj
   */
  void print_REMD(std::ostream &os, util::Replica_Data const & replica_data,
		  simulation::Parameter const & param,
		  int reeval = 0);
			    
   /** 
    * print RAMD information
    */
   void print_RAMD(std::ostream &os, configuration::Configuration const & conf, double lambda);

     /**
   * Print the SASA block.
   */
  void print_sasa(std::ostream &os, topology::Topology const & topo,
                    configuration::Configuration const & conf,
                    simulation::Simulation const & sim,
		    std::string const title);
  /**
   * Print sasa and volume averages and fluctuation.
   */
  void print_sasa_avg(std::ostream &os, std::vector<double> const &s,
                      std::vector<double> const &v, double &stot, double &vtot,
                      std::string const title, int volume);
  void print_sasa_fluct(std::ostream &os, std::vector<double> const &s,
                        std::vector<double> const &v, double &stot, double &vtot,
                        std::string const title, int volume);

  // delete after testing
  void print_forces(std::ostream &os, topology::Topology const & topo,
                    configuration::Configuration const & conf,
                    simulation::Simulation const & sim);

} // io

#endif
