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
 * @file perturbed_nonbonded_innerloop.h
 * eds-perturbed inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_INNERLOOP_H
#define INCLUDED_PERTURBED_NONBONDED_INNERLOOP_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Innerloop
   * eds-perturbed non bonded inner loop.
   */
  template<typename t_interaction_spec,
	   typename t_perturbation_details>
  class Perturbed_Nonbonded_Innerloop:
    public Eds_Nonbonded_Term
  {
  public:

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;
    
    /**
     * Constructor
     */
    explicit Perturbed_Nonbonded_Innerloop(Nonbonded_Parameter &nbp): m_param(&nbp){}
    
    /**
     * eds-perturbed interaction
     */
    void perturbed_lj_crf_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i, unsigned int j,
     Storage &storage,
     Periodicity_type const & periodicity,
    simulation::Simulation & sim
     );

    /**
     * eds-perturbed 1-4 interaction
     * (always shortrange)
     */
    void perturbed_one_four_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     unsigned int i, unsigned int j,
     Periodicity_type const & periodicity,
     simulation::Simulation & sim);
    
    /**
     * eds-perturbed RF interaction (solute).
     * (always shortrange)
     */
    void eds_RF_excluded_interaction_innerloop
    ( topology::Topology & topo,
      configuration::Configuration & conf,
      std::map<unsigned int, topology::EDS_Perturbed_Atom>::const_iterator const & mit,
      Periodicity_type const & periodicity, simulation::Simulation & sim);

     /**
     * perturbed RF interaction (solute).
     * (always shortrange)
     */
    void perturbed_RF_excluded_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
        std::map<unsigned int, topology::Perturbed_Atom>::const_iterator const & mit,
        Periodicity_type const & periodicity, simulation::Simulation & sim
    );

    /**
     * Common part of eds rf and pertubed atoms rf
     */
    void eds_perturbed_RF_exclusions_loop
    (topology::Topology & topo, 
    configuration::Configuration & conf, 
    int atom_i, const topology::excl_cont_t::value_type &exclusions, 
    Periodicity_type const & periodicity, simulation::Simulation & sim);

    /**
     * Calculation of the perturbed electric field (polarisation)
     */
    void perturbed_electric_field_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i, unsigned int j, math::Vec &e_eli, math::Vec &e_elj,
     Periodicity_type const & periodicity
    );
    
    /**
     * Calculation of the perturbed self energy (polarisation)
     */
    void perturbed_self_energy_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i,
     Storage & storage,
     Periodicity_type const & periodicity
    );

    
  protected:

    Nonbonded_Parameter * m_param;
    
  };
  
} // interaction

#include "perturbed_nonbonded_innerloop.cc"

#endif
