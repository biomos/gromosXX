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
 * @file cuda_pairlist_algorithm.h
 * CUDA accelerated pairlist algorithm
 */

#pragma once

#include "cuda_pairlist_algorithm_impl.h"

namespace math
{
  template<math::boundary_enum>
  class Periodicity;
}

namespace interaction
{
  class Pairlist;
  
  // template<typename t_interaction_spec>
  // class Nonbonded_Innerloop;
  
  // template<typename t_interaction_spec, typename t_perturbation_details>
  // class Perturbed_Nonbonded_Innerloop; 
  
  /**
   * @class CUDA_Pairlist_Algorithm
   * create an atomic pairlist with a
   * chargegroup based or atom based
   *  cut-off criterion.
   */

  template <typename Backend = util::cpuBackend>
  class CUDA_Pairlist_Algorithm : 
    public Pairlist_Algorithm, private algorithm::AlgorithmB<Backend>
  {
  public:
    /**
     * Constructor.
     */
    CUDA_Pairlist_Algorithm();

    /**
     * Destructor.
     */
    virtual ~CUDA_Pairlist_Algorithm() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      if (!quiet)
        os << "\tcuda pairlist algorithm\n";
      
      // initialize Cuda variables
      // maybe also copy simulation constants


  
      return 0;
    };

    /**
     * prepare the pairlists
     */    
    virtual int prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim);

    /**
     * update the pairlist
     */
    virtual void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			interaction::PairlistContainer & pairlist,
			unsigned int begin, unsigned int end,
			unsigned int stride);

    /**
     * update the pairlist, separating perturbed and non-perturbed interactions
     */
    virtual void update_perturbed(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
                                  interaction::PairlistContainer & pairlist,
				  interaction::PairlistContainer & perturbed_pairlist,
				  unsigned int begin, unsigned int end, 
				  unsigned int stride);

    bool excluded_solute_pair(topology::Topology & topo,
			      unsigned int i, unsigned int j);

    void set_cutoff(double const cutoff_short, double const cutoff_long)
    {
      m_cutoff_long = cutoff_long;
      m_cutoff_short = cutoff_short;
      m_cutoff_short_2 = cutoff_short * cutoff_short;
      m_cutoff_long_2  = cutoff_long * cutoff_long;
    }

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      os << "            "
	 << std::setw(32) << std::left << "solv - solv pairlist"
	 << std::setw(20) << m_solvent_solvent_timing << "\n";
    }
      
  private:
    CUDA_Pairlist_Algorithm_Impl<Backend> m_impl;
    /**
     * the chargegroup center of geometries. GPU allocated version
     */
    math::VArray m_cg_cog;
    /**
     * squared shortrange cutoff.
     */
    double m_cutoff_short_2;
    /**
     * squared longrange cutoff.
     */
    double m_cutoff_long_2;
    /**
     * longrange cutoff.
     */
    double m_cutoff_long;
    /**
     * shortrange cutoff.
     */
    double m_cutoff_short;
    /**
     * solvent - solvent pairlist 
     */
    double m_solvent_solvent_timing;
    
    math::Vec m_half_box;
    math::Vec m_box;

  };
} // interaction


