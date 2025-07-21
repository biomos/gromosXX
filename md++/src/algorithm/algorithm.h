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

#ifndef INCLUDED_ALGORITHM_H
#define INCLUDED_ALGORITHM_H

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
   * @class IAlgorithm
   * base class
   */
  class IAlgorithm
  {
  public:
    /**
     * Constructor.
     * @param name of the algorithm.
     */
    IAlgorithm(std::string name) : name(name), m_timer(name) {}

    /**
     * Destructor.
     */
    virtual ~IAlgorithm() {}
    
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
   * @class AlgorithmT
   * template class
   */
template <typename Backend = util::cpuBackend>
  class AlgorithmT : public IAlgorithm {
    /**
     * @brief compile time check for CUDA enabled
     * 
     */
    static_assert(util::is_valid_algorithm_backend<Backend>::value,
                "This backend is not supported in the current build configuration.");
  public:
    /**
     * Default constructor.
     */
    // Constructor forwarding name to base IAlgorithm
    explicit AlgorithmT(const std::string& name) : IAlgorithm(name) {}
    

    /**
     * Destructor.
     */
    virtual ~AlgorithmT() override {}
    
    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false) override = 0;
    // { return 0; }
    
    /**
     * apply the algorithm
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim) override {return 0;}
  };

  /**
   * @brief Allow use of Algorithm directly - defaults to AlgorithmT<util::cpuBackend>
   * 
   */
  using Algorithm = AlgorithmT<util::cpuBackend>;

  /**
   * @brief Create a backend-aware algorithm instance (GPU if available and supported, otherwise CPU)
   *
   * @tparam AlgT The algorithm template
   * @param os Output stream for logging
   * @return std::unique_ptr<IAlgorithm>
   */
//   template <template <typename> class AlgT>
//   std::unique_ptr<IAlgorithm> make_algorithm(std::ostream &os = std::cout) {
// #ifdef USE_CUDA
//     if constexpr (has_gpu_backend_v<AlgT>) {
//       if (gpu::CudaManager::is_enabled()) {
//         return std::make_unique<AlgT<util::gpuBackend>>(os);
//       }
//     }
// #endif
//     return std::make_unique<AlgT<util::cpuBackend>>(os);
//   }

  template <template <typename> class AlgT, typename... Args>
  IAlgorithm* make_algorithm(
                            simulation::Simulation & sim, 
                            Args&&... args) {
    if constexpr (util::has_gpu_backend_v<AlgT>) {
      if (sim.cuda_enabled()) {
        return new AlgT<util::gpuBackend>(std::forward<Args>(args)...);
      }
    }
    return new AlgT<util::cpuBackend>(std::forward<Args>(args)...);
  }
}
#endif

