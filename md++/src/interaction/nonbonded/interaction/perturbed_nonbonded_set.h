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
 * @file perturbed_nonbonded_set.h
 * the (eds-perturbed) non bonded interactions for a set of atoms:
 * Lennard-Jones and Coulomb interactions
 * eds-perturbed and non-eds-perturbed.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_SET_H
#define INCLUDED_PERTURBED_NONBONDED_SET_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Set
   * calculates the (eds-perturbed) nonbonded interactions.
   */
  class Perturbed_Nonbonded_Set : public Nonbonded_Set
  {
  public:    
    /**
     * Constructor.
     */
    Perturbed_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg,
			    Nonbonded_Parameter & param,
			    int rank, int num_threads);

    /**
     * Destructor
     */
    virtual ~Perturbed_Nonbonded_Set() {}
    
    /**
     * initialize some things
     */
    virtual int init(topology::Topology const & topo,
		     configuration::Configuration const & conf,
		     simulation::Simulation const & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * store results in configuration
     */
    virtual int update_configuration(topology::Topology const & topo,
				     configuration::Configuration & conf,
				     simulation::Simulation const & sim);

    /**
     * calculate the hessian for a given atom.
     */
    virtual int calculate_hessian(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  unsigned int atom_i, unsigned int atom_j,
				  math::Matrix & hessian);
   
    PairlistContainer & perturbed_pairlist() { return m_perturbed_pairlist; }
    PairlistContainer const & perturbed_pairlist()const { return m_perturbed_pairlist; }

  protected:

    PairlistContainer m_perturbed_pairlist;

    Perturbed_Nonbonded_Outerloop m_perturbed_outerloop;

    Perturbed_Nonbonded_Pair m_perturbed_pair;
    
  };
  
} // interaction

#endif
