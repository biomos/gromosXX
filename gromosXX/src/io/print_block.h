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
  inline std::ostream &
  print_MULTIBATH_COUPLING(std::ostream &os, simulation::Multibath const &bath);

  /** 
   * Print the MULTIBATH block.
   */
  inline std::ostream & 
  print_MULTIBATH(std::ostream &os, simulation::Multibath const &bath,
		  simulation::Energy const &energy);

  /** 
   * Print the DEGREES OF FREEDOM block.
   */
  inline std::ostream & 
  print_DEGREESOFFREEDOM(std::ostream &os, simulation::Multibath const &bath);

  /**
   * Print the PCOUPLE block.
   */
  inline std::ostream &
  print_PCOUPLE(std::ostream &os,
		bool calc, int ntp, math::Matrix pres0, double comp, 
		double tau, interaction::virial_enum vir);

  /**
   * Print the PRESSURE block.
   */
  template<math::boundary_enum b>
  inline std::ostream &
  print_PRESSURE(std::ostream &os, simulation::System<b> const & sys);
  
  /**
   * Print the ENERGY block.
   */
  inline std::ostream &
  print_ENERGY(std::ostream &os, simulation::Energy const &e,
	       std::vector<size_t> const & energy_groups,
	       std::string const title = "ENERGIES");
  
  /**
   * Print the TIMESTEP block.
   */
  inline std::ostream &
  print_TIMESTEP(std::ostream &os, double const steps, double const time);

  /**
   * Print the CENTREOFMASS block.
   */
  inline std::ostream &
  print_CENTREOFMASS(std::ostream &os, double const ekin_trans, double const ekin_rot);
  
} // io

#include "print_block.tcc"

#endif
