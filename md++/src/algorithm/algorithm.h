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
#include "gpu/mirror_fields.h"
#include "io/message.h"

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
     * @brief Which Configuration-mirror fields (gpu::MirrorField bits)
     * this algorithm's apply() might have written directly on the CPU
     * side, for CudaManager's GPU-mirror freshness tracking.
     * Conservative default: assume everything might have changed, so
     * Algorithm_Sequence::run() always invalidates the GPU mirror after
     * algorithms that don't override this -- correct by default for
     * every existing (CPU-side) algorithm with zero changes needed.
     * Only the few algorithms that manage their own GPU-mirror
     * freshness through sim.cuda() (currently just
     * Leap_Frog_Velocity<gpuBackend>) need to narrow this, so the
     * sequence's default invalidation doesn't immediately erase what
     * they just marked fresh.
     */
    virtual unsigned gpu_mirror_touches() const { return gpu::MIRROR_ALL; }

    /**
     * @brief Does this algorithm's apply() need any *other* algorithm's
     * GPU-deferred work (see finalize_gpu_step() below) resolved on the
     * CPU before it runs? Default false. Algorithm_Sequence::run()
     * checks this before each apply() and flushes any pending deferred
     * work first -- e.g. Energy_Calculation reads conf.old().energies.
     * kinetic_energy (calculate_totals()), which Temperature_
     * Calculation<gpuBackend> may have only queued asynchronously on
     * the GPU rather than finished computing on the host yet.
     */
    virtual bool needs_finalized_gpu_state() const { return false; }

    /**
     * @brief Did this algorithm's most recent apply() leave GPU work
     * queued (kernels launched, no CPU sync) whose result some later
     * algorithm this same step might need? Default false -- the
     * overwhelming majority of algorithms either don't touch the GPU at
     * all or already sync within their own apply(). Only an algorithm
     * that overrides this to possibly return true needs finalize_gpu_
     * step() to do real work; Algorithm_Sequence::run() uses this to
     * decide which algorithms to add to (and later flush from) its
     * per-step pending list, without needing to call finalize_gpu_step()
     * unconditionally on everything every step.
     */
    virtual bool has_pending_gpu_finalize() const { return false; }

    /**
     * @brief Resolve GPU work this algorithm's apply() left queued
     * (sync + whatever host-side bookkeeping depends on the result).
     * Called by Algorithm_Sequence::run(): once before any later
     * algorithm in the same step whose needs_finalized_gpu_state() is
     * true, and once more, unconditionally, as a safety net at the very
     * end of run() (so nothing is ever silently left unresolved just
     * because no consumer happened to run this step). Default no-op.
     */
    virtual void finalize_gpu_step(topology::Topology &,
                                    configuration::Configuration &,
                                    simulation::Simulation &) {}

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

  namespace detail {
    /**
     * @brief Mixed CPU/GPU diagnostic (see Interaction::is_gpu_native()'s
     * doc comment for the Interaction-side equivalent): a single warning,
     * emitted once at construction, when GPU acceleration is active but
     * this particular Algorithm has no GPU backend at all and is running
     * on CPU instead. Shared by every make_algorithm/make_unique_algorithm
     * overload below rather than duplicated at each one.
     */
    inline void warn_if_cpu_fallback(const simulation::Simulation & sim,
                                      const std::string & name) {
      if (sim.param().gpu.accelerator == simulation::gpu_cuda) {
        io::messages.add(
            "Algorithm '" + name + "' has no GPU implementation and will "
            "run on CPU -- if it reads or writes positions/velocities/"
            "force, this adds extra CPU<->GPU synchronization every step.",
            "make_algorithm", io::message::warning);
      }
    }
  }

  /**
   * @brief Create a backend-aware algorithm instance (GPU if available and supported, otherwise CPU)
   *
   * @tparam AlgT The algorithm template
   * @param Args... Arguments to be passed to the constructor
   * @return Algorithm*
   */
  template <template <typename> class AlgT, typename... Args>
  Algorithm* make_algorithm(const simulation::Simulation & sim,
                            Args&&... args) {
    static_assert(std::is_base_of_v<Algorithm, AlgT<util::cpuBackend>>,
                  "AlgT must derive from Algorithm");
    if constexpr (util::has_gpu_backend_v<AlgT>) {
      if (sim.param().gpu.accelerator == simulation::gpu_cuda) {
        return new AlgT<util::gpuBackend>(std::forward<Args>(args)...);
      }
    }
    Algorithm * alg = new AlgT<util::cpuBackend>(std::forward<Args>(args)...);
    detail::warn_if_cpu_fallback(sim, alg->name);
    return alg;
  }

  /**
   * @brief Legacy version of make_algorithm for non-templated classes --
   * these have no Backend parameter at all, so under GPU acceleration
   * they unconditionally warn (see detail::warn_if_cpu_fallback()).
   *
   * @tparam Alg The algorithm class
   * @param sim Simulation object
   * @param args Arguments to be passed to the constructor
   * @return Alg*
   */
  template <class Alg, typename... Args>
  Alg* make_algorithm(const simulation::Simulation & sim,
                      Args&&... args) {
    Alg * alg = new Alg(std::forward<Args>(args)...);
    detail::warn_if_cpu_fallback(sim, alg->name);
    return alg;
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
    std::unique_ptr<Algorithm> alg = std::make_unique<AlgT<util::cpuBackend>>(std::forward<Args>(args)...);
    detail::warn_if_cpu_fallback(sim, alg->name);
    return alg;
  }
}
