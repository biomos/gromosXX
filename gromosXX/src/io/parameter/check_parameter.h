/**
 * @file check_parameter.h
 * check parameters
 */

#ifndef INCLUDED_CHECK_PARAMETER_H
#define INCLUDED_CHECK_PARAMETER_H

namespace io
{
  /**
   * cross checks on the parameters
   */
  int check_parameter(simulation::Simulation &sim,topology::Topology const & topo);
  
  /**
   * does basic cross checks on parameters from different blocks
   * cross-checks of parameters within one block should be done in the 
   * read_XXX functions in @ref In_Parameter
   */
  int simple_crosschecks(simulation::Simulation &sim);
  
  /**
   * does basic cross checks on the parameters
   * checks each feature against all others
   * every allowed combination has to be unlocked
   */
  int check_features(simulation::Simulation &sim,topology::Topology const &topo);
}

#endif
