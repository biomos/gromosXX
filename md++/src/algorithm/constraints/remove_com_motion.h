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
 * @file remove_com_motion.h
 * remove center of mass translational and angular momentum
 */

#pragma once

#include "gpu/cuda/manager/cuda_manager.h"

namespace algorithm
{
  /**
   * @class Remove_COM_Motion
   * implements centre of mass motion removal
   * input switches determine whether to remove
   * translational or angular centre of mass
   * momentum (or both).
   * For periodic boundary conditions, only
   * translational centre of mass motion is
   * removed.
   * It is recommended to remove it at every step.
   *
   * @todo check if there is a Gromos96 bug in cenmas.
   * seems like absolute positions instead of
   * relative to centre of mass are taken in
   * angular momentum calculation.
   */
  template<typename Backend = util::cpuBackend>
  class Remove_COM_Motion : public Algorithm, private AlgorithmB<Backend>
  {
  public:
    /**
     * Specify supported backends
     * 
     */
    template <typename B>
    static constexpr bool is_supported_backend =
            std::is_same_v<B, util::cpuBackend>
        ||  std::is_same_v<B, util::gpuBackend>
      ;
    /**
     * Constructor.
     */
    Remove_COM_Motion(std::ostream & os = std::cout) : Algorithm("RemoveCOMMotion"), os(os) {}

    /**
     * Destructor.
     */
    virtual ~Remove_COM_Motion() {}
    
    /**
     * apply COM removal.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

    /**
     * calculate and remove translational centre of mass motion
     */
    double remove_com_translation(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  bool remove_trans = true);
    
    /**
     * calculate and remove angular centre of mass motion
     */
    double remove_com_rotation(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       bool remove_rot = true);

    /**
     * add centre of mass rotation
     */
    double add_com_rotation(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    math::Vec com_L);

    // gpuBackend: remove_com_translation()/remove_com_rotation() manage
    // POS/VEL freshness themselves (configuration_view()/mark_gpu_dirty()
    // in remove_com_motion_gpu.cc) -- narrowing to 0 stops Algorithm_
    // Sequence::run()'s default before/after hooks from forcing an
    // eager publish on every step (this is the FIRST algorithm in
    // create_md_sequence.cc's sequence; leaving it at the base class
    // default MIRROR_ALL forced a full mirror round-trip every single
    // step regardless of comtransrot's skip-step cadence or of what any
    // downstream deferred-sync algorithm had managed to avoid -- see
    // PERFORMANCE.md's "Architecture direction" section). On the
    // (common) steps where apply() is a no-op (skip_step not due),
    // this mask means the framework does nothing around it at all, not
    // even an unnecessary no-op flush check. cpuBackend: matches the
    // base class default (MIRROR_ALL).
    virtual unsigned gpu_mirror_touches() const override {
      if constexpr (std::is_same_v<Backend, util::gpuBackend>)
        return 0u;
      else
        return gpu::MIRROR_ALL;
    }

  protected:
    std::ostream & os;
  };

  // extern template class Remove_COM_Motion<util::cpuBackend>;
  // #ifdef USE_CUDA
  // extern template class Remove_COM_Motion<util::gpuBackend>;
  // #endif
} //algorithm

// #include "remove_com_motion_cpu.tcc"
// #ifdef USE_CUDA
//   #include "remove_com_motion_gpu.tcc"
// #endif
