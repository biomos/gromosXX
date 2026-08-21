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
 * @file lattice_shift.h
 * keeping track of lattice shifts
 */

#ifndef INCLUDED_LATTICE_SHIFT_H
#define INCLUDED_LATTICE_SHIFT_H

#include <type_traits>
namespace algorithm
{
  /**
   * @class Lattice_Shift_Tracker
   * keeps track of lattice shifts
   */
  template <typename Backend = util::cpuBackend>
  class Lattice_Shift_Tracker : public Algorithm, private AlgorithmB<Backend>
  {
  public:
    template <typename B>
    static constexpr bool is_supported_backend =
        std::is_same_v<B, util::cpuBackend> ||
        std::is_same_v<B, util::gpuBackend>;

    /**
     * Constructor.
     */
    Lattice_Shift_Tracker() : Algorithm("Lattice_Shift_Tracker") {}

    /**
     * Destructor.
     */
    virtual ~Lattice_Shift_Tracker(){}

    /**
     * put CG into box and keep track of shift
     */
    virtual int apply(topology::Topology &topo,
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

    // gpuBackend: manages POS and LATTICE_SHIFT freshness itself
    // (configuration_view()/mark_gpu_dirty() in lattice_shift_gpu.cc),
    // same reasoning as Leap_Frog_Position<gpuBackend> -- narrowing to
    // 0 stops Algorithm_Sequence::run()'s default before/after hooks
    // from forcing an eager publish or erasing the freshness apply()
    // just set. This is the actual point of the GPU port: RemoveCOM
    // Motion/Lattice_Shift_Tracker together used to be the two
    // default-MIRROR_ALL algorithms sitting at the very top of every
    // step, forcing a full mirror round-trip regardless of what any
    // downstream deferred-sync algorithm (Leap_Frog_Position, Berendsen
    // _Barostat, ...) had managed to avoid -- see PERFORMANCE.md.
    // cpuBackend: matches the base class default (MIRROR_ALL).
    virtual unsigned gpu_mirror_touches() const override {
      if constexpr (std::is_same_v<Backend, util::gpuBackend>)
        return 0u;
      else
        return gpu::MIRROR_ALL;
    }

  protected:
    template<math::boundary_enum b>
    void _apply(topology::Topology &topo,
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);
  };
}
#endif	/* INCLUDED_LATTICE_SHIFT_H */

