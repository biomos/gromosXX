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
  /**
   * @class CUDA_Pairlist_Algorithm
   * create an atomic pairlist on GPU
   * with a chargegroup based or atom
   * based cut-off criterion.
   */

  template <typename Backend = util::cpuBackend>
  class CUDA_Pairlist_Algorithm : public Pairlist_Algorithm, private algorithm::AlgorithmB<Backend>
    /**
     * We have to inherit from Pairlist_Algorithm to have a common class pointer
     */
    
  {
  public:
    /**
     * Constructor.
     */
    CUDA_Pairlist_Algorithm();

    /**
     * destructor.
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
     * prepare the pairlist(s).
     */
    virtual int prepare(topology::Topology & topo,
                        configuration::Configuration & conf,
                        simulation::Simulation &sim);

    /**
     * @brief 
     * 
     * @tparam PairlistContainerType can be interaction::PairlistContainer or gpu::PairlistContainer
     * @param topo 
     * @param conf 
     * @param sim 
     * @param pairlist 
     * @param begin 
     * @param end 
     * @param stride 
     */
    // using PairlistContainerType = interaction::PairlistContainer;
    // template <typename PairlistContainerType>
    virtual void update(topology::Topology & topo,
                        configuration::Configuration & conf,
                        simulation::Simulation &sim,
                        interaction::PairlistContainer &pairlist,
                        unsigned int begin, unsigned int end, 
                        unsigned int stride);

    virtual void update_perturbed(topology::Topology & topo,
                                  configuration::Configuration & conf,
                                  simulation::Simulation & sim,
                                  interaction::PairlistContainer & pairlist,
                                  interaction::PairlistContainer & perturbed_pairlist,
                                  unsigned int begin, unsigned int end, 
                                  unsigned int stride) {
      std::cerr <<__FILE__ << ":" << __LINE__ <<
      " do not use this overload of update_perturbed on a CUDA pairlist algorithm" << std::endl;
      // assert(false);
    };

  private:
    CUDA_Pairlist_Algorithm_Impl<Backend> m_impl;
  };
} // interaction


