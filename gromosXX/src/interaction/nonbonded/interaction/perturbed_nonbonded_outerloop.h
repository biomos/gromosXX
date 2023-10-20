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
 * @file perturbed_nonbonded_outerloop.h
 * the eds-perturbed non bonded outer loop:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_OUTERLOOP_H
#define INCLUDED_PERTURBED_NONBONDED_OUTERLOOP_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Outerloop
   * outerloop for the eds-perturbed nonbonded interactions.
   */
  class Perturbed_Nonbonded_Outerloop
  {
  public:    

    /**
     * Constructor.
     * @param nbp the nonbonded parameters
     */
    Perturbed_Nonbonded_Outerloop(Nonbonded_Parameter & nbp);

    /**
     * calculate the eds-perturbed interactions.
     */
    void perturbed_lj_crf_outerloop(topology::Topology & topo,
				    configuration::Configuration & conf,
				    simulation::Simulation & sim,
				    Pairlist const & pairlist,
				    Storage & storage);  

    /**
     * calculate the eds-perturbed 1,4-interactions.
     */
    void perturbed_one_four_outerloop(topology::Topology & topo,
				      configuration::Configuration & conf,
				      simulation::Simulation & sim,
				      Storage & storage);

    /**
     * calculate the eds-perturbed RF contributions for excluded atoms.
     */
    void eds_RF_excluded_outerloop(topology::Topology & topo,
					 configuration::Configuration & conf,
					 simulation::Simulation & sim,
					 Storage & storage);

    /**
     * calculate the perturbed self energy (polarisation).
     */
    void perturbed_self_energy_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
			    Storage & storage);

    /**
     * calculate the perturbed electric field (polarisation).
     */
    void perturbed_electric_field_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
		            PairlistContainer const & pairlist,
                            PairlistContainer const & perturbed_pairlist,
                            Storage & storage,
                            Storage & storage_lr,
                            int rank);

  private:
    Nonbonded_Parameter & m_param;

    /**
     * calculate the eds-perturbed interactions.
     */
    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _perturbed_lj_crf_outerloop(topology::Topology & topo,
				     configuration::Configuration & conf,
				     simulation::Simulation & sim,
				     Pairlist const & pairlist,
				     Storage & storage);
       
    /**
     * calculate the eds-perturbed 1,4-interactions.
     */
    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _perturbed_one_four_outerloop(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim,
				       Storage & storage);

    /**
     * calculate the eds-perturbed RF contributions for excluded atoms.
     */
    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _eds_RF_excluded_outerloop(topology::Topology & topo,
					  configuration::Configuration & conf,
					  simulation::Simulation & sim,
					  Storage & storage);

    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _perturbed_self_energy_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
			    Storage & storage);

    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _perturbed_electric_field_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
		            PairlistContainer const & pairlist,
                            PairlistContainer const & perturbed_pairlist,
			    Storage & storage,
                            Storage & storage_lr,
                            int rank);
  };
  
} // interaction

#endif
