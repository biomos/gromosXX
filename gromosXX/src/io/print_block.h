/**
 * @file print_block.h
 * routines to print out the various blocks.
 */

#ifndef INCLUDED_PRINT_BLOCK_H
#define INCLUDED_PRINT_BLOCK_H

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
		       configuration::Energy const &energy);

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
		     double tau, math::virial_enum vir);

  /**
   * Print the PRESSURE block.
   */
  void print_PRESSURE(std::ostream &os, configuration::Configuration const & conf);
  
  /**
   * Print the ENERGY block.
   */
  void print_ENERGY(std::ostream &os, configuration::Energy const &e,
		    std::vector<size_t> const & energy_groups,
		    std::string const title = "ENERGIES");
  
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
  
} // io

#endif
