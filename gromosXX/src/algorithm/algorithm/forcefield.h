/**
 * @file forcefield.h
 * global functions to create a forcefield.
 */

#ifndef INCLUDED_GLOBAL_FORCEFIELD_H
#define INCLUDED_GLOBAL_FORCEFIELD_H

namespace algorithm
{

  /**
   * create a Gromos96 like forcefield.
   */
  template<typename t_simulation, typename t_interaction_spec>
  void G96_Forcefield(interaction::Forcefield<
		      t_simulation, t_interaction_spec> &ff,
		      t_simulation & sim,
		      io::InTopology &topo,
		      io::InInput &input,
		      io::Argument &args);

  /**
   * create a perturbed Gromos96 like forcefield.
   */
  template<typename t_simulation, typename t_interaction_spec>
  void Perturbed_G96_Forcefield(interaction::Forcefield<
				t_simulation, t_interaction_spec> &ff,
				t_simulation & sim,
				io::InTopology &topo,
				io::InInput &input,
				io::Argument &args);
  
} // algorithm

#include "forcefield.tcc"

#endif
