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
 * @file berendsen_barostat.h
 * berendsen barostat
 *
 * Templated on Backend:
 *   - util::cpuBackend : existing per-atom loop, all four pcouple.scale
 *     modes (isotropic/anisotropic/full_anisotropic/semi_anisotropic).
 *   - util::gpuBackend : the O(9) mu-matrix computation and box scaling
 *     stay on the host (trivial, box is a 3x3 struct); the O(num_atoms)
 *     position scaling runs as a single GPU kernel, GPU-resident (no
 *     eager CPU publish -- see gpu_mirror_touches() below). Only
 *     isotropic/anisotropic/semi_anisotropic are supported on GPU --
 *     full_anisotropic is refused with an error (see berendsen_
 *     barostat_gpu.cc's apply() for why: the CPU reference has a
 *     pre-existing fallthrough bug for this mode that isn't this
 *     branch's to silently replicate or silently fix).
 */

#pragma once

namespace algorithm
{

  /**
   * @class Berendsen_Barostat
   * the Berendsen barostat.
   */
  template <typename Backend = util::cpuBackend>
  class Berendsen_Barostat : public Algorithm, private AlgorithmB<Backend>
  {
  public:
    template <typename B>
    static constexpr bool is_supported_backend =
        std::is_same_v<B, util::cpuBackend> ||
        std::is_same_v<B, util::gpuBackend>;

    /**
     * Constructor.
     */
    Berendsen_Barostat() : Algorithm("BerendsenBarostat") {}
    /**
     * Destructor.
     */
    virtual ~Berendsen_Barostat() {}

    /**
     * apply weak coupling.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    // gpuBackend: box scaling still writes conf.current().box directly
    // on the host every call (cheap, O(9)) -- BOX stays in the mask so
    // the after-invalidate lets a later GPU reader (e.g. the pairlist)
    // resync it. POS is deliberately left out: apply() manages its
    // freshness itself via mark_gpu_dirty() (GPU-resident scaling, no
    // CPU round trip), and including it here would make Algorithm_
    // Sequence::run()'s default before-apply() flush_gpu_dirty() defeat
    // that by eagerly publishing it right before we scale it anyway.
    // cpuBackend: matches the base class default (MIRROR_ALL) --
    // nothing GPU-resident to protect, box+pos are plain CPU writes.
    virtual unsigned gpu_mirror_touches() const override {
      if constexpr (std::is_same_v<Backend, util::gpuBackend>)
        return gpu::MIRROR_BOX;
      else
        return gpu::MIRROR_ALL;
    }

    /**
     * init
     */
    virtual int init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false)
    {
      // os << "Berendsen barostat\n";
      return 0;
    };

  private:

  };

} // algorithm
