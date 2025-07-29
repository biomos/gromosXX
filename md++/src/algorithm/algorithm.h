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
 * @file algorithm.h
 * base class for algorithms
 */

#pragma once

#include "simulation/simulation.h"

namespace configuration
{
  class Configuration;
}
namespace topology
{
  class Topology;
}
namespace simulation
{
  class Simulation;
}
namespace util
{
  class Algorithm_Timer;
}

namespace algorithm
{
  /**
   * @class Algorithm
   * base class
   */
  class Algorithm
  {
  public:
    /**
     * Constructor.
     * @param name of the algorithm.
     */
    Algorithm(std::string name) : name(name), m_timer(name) {}

    /**
     * Destructor.
     */
    virtual ~Algorithm() {}
    
    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false) = 0;
    // { return 0; }
    
    /**
     * apply the algorithm
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim) {return 0;}

    /**
     * name of the algorithm
     */
    std::string name;

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      m_timer.print(os);
    }
    /**
     * const accessor to timer
     */
    const util::Algorithm_Timer & timer() const {
      return m_timer;
    }
    /**
     * accessor to timer
     */
    util::Algorithm_Timer & timer() {
      return m_timer;
    }
    /**
     * accessor to timer
     */
    void timer(util::Algorithm_Timer &t) {
      m_timer = t;
    }
    
  protected:
    /**
     * store time used in algorithm.
     */
    util::Algorithm_Timer m_timer;
  };

  /**
   * @struct AlgorithmB
   * @brief Base class template for backend validation
   * 
   */
  template <typename Backend, typename Enable = void>
  struct AlgorithmB;

  /**
   * @struct AlgorithmB
   * @details Compile-time helper struct to check for supported backends
   */
  template <typename Backend>
  struct AlgorithmB<
      Backend,
      std::enable_if_t<
           std::is_same_v<Backend, util::cpuBackend>
        || std::is_same_v<Backend, util::gpuBackend>>
  > {
    /**
     * @brief Specify supported backends. This way we control 
     * make_algorithm factories. Supported are
     * util::cpuBackend or util::gpuBackend
     * 
     */
    template <typename B>
    static constexpr bool is_supported_backend =
            std::is_same_v<B, util::cpuBackend>
          /* uncomment this line in your algorithm to allow gpu backend for the algorithm */
        // ||  std::is_same_v<B, util::gpuBackend>
      ;
  };

  /**
   * @brief Create a backend-aware algorithm instance (GPU if available and supported, otherwise CPU)
   *
   * @tparam AlgT The algorithm template
   * @param Args... Arguments to be passed to the constructor
   * @return Algorithm*
   */
  template <template <typename> class AlgT, typename... Args>
  Algorithm* make_algorithm(simulation::Simulation & sim, 
                            Args&&... args) {
    static_assert(std::is_base_of_v<Algorithm, AlgT<util::cpuBackend>>,
                  "AlgT must derive from Algorithm");
    if constexpr (util::has_gpu_backend_v<AlgT>) {
      if (sim.param().gpu.accelerator == simulation::gpu_cuda) {
        return new AlgT<util::gpuBackend>(std::forward<Args>(args)...);
      }
    }
    return new AlgT<util::cpuBackend>(std::forward<Args>(args)...);
  }

  /**
   * @brief Legacy version of make_algorithm for non-templated classes
   * 
   * @tparam Alg The algorithm class
   * @param sim Simulation object
   * @param args Arguments to be passed to the constructor
   * @return Alg* 
   */
  template <class Alg, typename... Args>
  Alg* make_algorithm(simulation::Simulation & sim, 
                      Args&&... args) {
    return new Alg(std::forward<Args>(args)...);
  }


  /**
   * @brief Create a backend-aware algorithm instance (GPU if available and supported, otherwise CPU)
   *
   * @tparam AlgT The algorithm template
   * @param Args... Arguments to be passed to the constructor
   * @return std::unique_ptr<Algorithm>
   */
  template <template <typename> class AlgT, typename... Args>
  std::unique_ptr<Algorithm> make_unique_algorithm(
                            simulation::Simulation & sim, 
                            Args&&... args) {
    if constexpr (util::has_gpu_backend_v<AlgT>) {
      if (sim.param().gpu.accelerator == simulation::gpu_cuda) {
        return std::make_unique<AlgT<util::gpuBackend>>(std::forward<Args>(args)...);
      }
    }
    return std::make_unique<AlgT<util::cpuBackend>>(std::forward<Args>(args)...);
  }
}
