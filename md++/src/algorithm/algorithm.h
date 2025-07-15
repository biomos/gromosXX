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
  template<typename Backend = util::cpuBackend>
  class AlgorithmT : public IAlgorithm, private util::BackendData<Backend> {
  public:
    // Constructors:
    template<typename B = Backend, typename = std::enable_if_t<std::is_same_v<B, util::cpuBackend>>>
    AlgorithmT(const std::string& name) : IAlgorithm(name) {}

    template<typename B = Backend, typename = std::enable_if_t<std::is_same_v<B, util::gpuBackend>>>
    AlgorithmT(std::shared_ptr<gpu::CudaManager> mgr, const std::string& name)
      : IAlgorithm(name) {
      this->cuda_manager_ = std::move(mgr);  // 'this->' needed due to dependent base
    }

#ifndef USE_CUDA
  static_assert(std::is_same_v<Backend, util::cpuBackend>,
                "This algorithm is allowed only with `util::cpu` unless CUDA is enabled (USE_CUDA)");
#endif

    // Accessor for GPU manager if needed
    auto& cuda_manager() {
      static_assert(std::is_same_v<Backend, util::gpuBackend>, "cuda_manager only valid for gpu backend");
      return this->cuda_manager_;
    }

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
}

#endif

