/**
 * @file print_block.h
 * routines to print out the various blocks.
 */

#ifndef INCLUDED_PRINT_BLOCK_H
#define INCLUDED_PRINT_BLOCK_H

namespace io
{
  /** 
   * Print the MULTIBATH block.
   */
  inline std::ostream & 
  print_MULTIBATH(std::ostream &os,simulation::Multibath const &bath);
  
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

} // io

#include "print_block.tcc"

#endif
