/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file qmmm_nonbonded_outerloop.h
 * the QMMM nonbonded outerloops.
 */

#ifndef INCLUDED_QMMM_NONBONDED_OUTERLOOP_H
#define INCLUDED_QMMM_NONBONDED_OUTERLOOP_H

#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"

namespace interaction
{
  
  /**
   * @class QMMM_Nonbonded_Outerloop
   * loops the nonbonded interactions...
   */
  class QMMM_Nonbonded_Outerloop : private Nonbonded_Outerloop
  {
  public:    
    /**
     * Constructor.
     */
    QMMM_Nonbonded_Outerloop(Nonbonded_Parameter & nbp);
    
    /**
     * calculate only lj interactions.
     */
    void lj_outerloop(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      Pairlist const & pairlist_solute,
                      Pairlist const & pairlist_solvent,
		      Storage & storage,
                      bool longrange,
                      util::Algorithm_Timer & timer,
                      bool master);

    /**
     * calculate the 1,4-interactions.
     */
    void one_four_outerloop(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    Storage & storage,
                            int rank, int size);

    /**
     * calculate the LJ exceptions
     */
    void lj_exception_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage,
            int rank, int size);
    
    /**
     * calculate the interaction for a given atom pair.
     * SLOW! as it has to create the periodicity...
     */
    int calculate_interaction(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim,
			      unsigned int atom_i, unsigned int atom_j,
			      math::Vec & force, 
			      double &e_lj, double &e_crf);

    /**
     * calculate the hessian for a given pair
     */
    int calculate_hessian(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  unsigned int atom_i, unsigned int atom_j,
			  math::Matrix & hessian,
			  PairlistContainer const & pairlist);

  private:
    template<typename t_interaction_spec>
    void _lj_outerloop(topology::Topology & topo,
                  configuration::Configuration & conf,
                  simulation::Simulation & sim,
                  Pairlist const & pairlist_solute,
                  Pairlist const & pairlist_solvent,
                  Storage & storage,
                  bool longrange, util::Algorithm_Timer & timer, bool master);

    template<typename t_interaction_spec>
    void _one_four_outerloop(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     Storage & storage,
                             int rank, int size);

    template<typename t_interaction_spec>
    void _lj_exception_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage,
            int rank, int size);
    
    template<typename t_interaction_spec>
    int _calculate_interaction(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       unsigned int atom_i, unsigned int atom_j,
			       math::Vec & force,
			       double & e_lj, double & e_crf);

    template<typename t_interaction_spec>
    int _calculate_hessian(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   unsigned int atom_i, unsigned int atom_j,
			   math::Matrix & hessian,
			   PairlistContainer const & pairlist);
  };
} // interaction

#endif
