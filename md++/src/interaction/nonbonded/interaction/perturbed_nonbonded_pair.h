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
 * @file perturbed_nonbonded_pair.h
 * perturbed inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_PAIR_H
#define INCLUDED_PERTURBED_NONBONDED_PAIR_H

namespace math
{
  template<math::boundary_enum b>
  class Periodicity;
}

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Pair
   * perturbed non bonded pairs.
   */
  class Perturbed_Nonbonded_Pair
  {
  public:
    
    /**
     * Constructor
     */
    Perturbed_Nonbonded_Pair(Nonbonded_Parameter &nbp);

    /**
     * calculate the perturbed pair contributions.
     */
    void perturbed_pair_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage);

  private:
    template<typename t_perturbation_details>
    void _perturbed_pair_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage);

    template<typename t_interaction_spec, typename t_perturbation_details>
    void _split_perturbed_pair_outerloop(topology::Topology & topo,
					 configuration::Configuration & conf,
					 simulation::Simulation & sim,
					 Storage & storage);

    /**
     * perturbed pairs.
     * (always shortrange)
     * NO RANGE FILTER FOR PERTURBED PAIRS ??
     * NO SCALING for PERTURBED PAIRS ??
     * NO MOLECULAR VIRIAL CONTRIBUTION ??
     */
    template<typename t_interaction_spec, typename t_perturbation_details, math::boundary_enum t_boundary_type>
    void perturbed_pair_interaction_innerloop
    ( topology::Topology & topo, configuration::Configuration & conf,
      simulation::Simulation & sim,
      std::vector<topology::perturbed_two_body_term_struct>
      ::const_iterator const &it,
      math::Periodicity<t_boundary_type> const & periodicity);
 
    Nonbonded_Parameter *m_param;
    Nonbonded_Term m_nonbonded_term;
    Eds_Nonbonded_Term m_perturbed_nonbonded_term;

  };
  
} // interaction

#endif
