/**
 * @file print_block.h
 * routines to print out the various blocks.
 */

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
  template<typename t_simulation>
  inline std::ostream &
  print_ENERGY(std::ostream &os, t_simulation const &sim);
  
} // io

#include "print_block.tcc"


  
