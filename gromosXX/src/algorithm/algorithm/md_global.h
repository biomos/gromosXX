/**
 * @file md_global.h
 * global functions to get an md simulation started.
 */

#ifndef INCLUDED_MD_GLOBAL_H
#define INCLUDED_MD_GLOBAL_H

namespace algorithm
{
  /**
   * perform an MD simulation.
   */
  int do_md(io::Argument &args);

  /**
   * helper function to do the md.
   */
  template<bool do_grid>
  int do_md_grid(io::Argument &args, io::InInput &input);

  /**
   * helper function to do the md.
   */
  template<bool do_grid, bool do_scaled>
  int do_md_scaled(io::Argument &args, io::InInput &input);

  /**
   * helper function to do the md.
   */
  template<bool do_grid, bool do_scaled, interaction::virial_enum do_virial>
  int do_md_virial(io::Argument &args, io::InInput &input);
  

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

#include "md_global.tcc"

#endif
