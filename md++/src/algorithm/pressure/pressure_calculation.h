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
 * @file pressure_calculation.h
 * pressure calculation
 *
 * Templated on Backend like Berendsen_Barostat: the O(9) matrix combine
 * (pressure_tensor = tot_cg_factor * (kinetic_energy_tensor +
 * 0.5*virial_tensor) * (2/volume)) is trivially cheap regardless of
 * backend and stays on the host either way -- what differs is how
 * virial_tensor is read. cpuBackend reads conf.old().virial_tensor
 * directly (always CPU-resident there). gpuBackend reads it through
 * configuration_view() (the same wait-on-fresh-producers protocol every
 * other GPU-native reader uses) instead of the generic flush_gpu_dirty()
 * path, which used to bundle this together with MIRROR_CONSTRAINT_FORCE
 * and could read a virial_tensor that a same-step CPU-side correction
 * (Molecular_Virial_Interaction) had written straight past the mirror --
 * see CUDA_Molecular_Virial_Interaction, which now keeps the mirror
 * itself authoritative so this direct read is correct.
 */

#ifndef INCLUDED_PRESSURE_CALCULATION_H
#define INCLUDED_PRESSURE_CALCULATION_H

namespace algorithm
{

  /**
   * @class Pressure_Calculation
   * pressure calculation.
   */
  template <typename Backend = util::cpuBackend>
  class Pressure_Calculation : public Algorithm, private AlgorithmB<Backend>
  {
  public:
    template <typename B>
    static constexpr bool is_supported_backend =
        std::is_same_v<B, util::cpuBackend> ||
        std::is_same_v<B, util::gpuBackend>;

    /**
     * Constructor.
     */
    Pressure_Calculation() : Algorithm("PressureCalculation") {}

    /**
     * Destructor.
     */
    virtual ~Pressure_Calculation() {}

    /**
     * apply the pressure calculation
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    /**
     * Reads conf.old().virial_tensor/kinetic_energy_tensor and writes
     * pressure_tensor -- all plain CPU-resident matrices computed by
     * this class or Molecular_Virial_Interaction/CUDA_Molecular_Virial_
     * Interaction, never GPU-mirror pos/vel/force/box data directly.
     * gpuBackend narrows away from the base class default (MIRROR_ALL):
     * it reads virial_tensor itself via configuration_view() (see this
     * class's own doc comment above), so the framework's default
     * before-apply() flush would be redundant work -- exactly the kind
     * of push this class doesn't need, muting the benefit of every
     * deferred-sync fix upstream of it (see PERFORMANCE.md's
     * "Architecture direction" section). cpuBackend matches the base
     * class default (MIRROR_ALL) -- virial_tensor is plain CPU data
     * there, no GPU mirror to protect against.
     */
    virtual unsigned gpu_mirror_touches() const override {
      if constexpr (std::is_same_v<Backend, util::gpuBackend>)
        return 0u;
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
      // os << "Pressure calculation\n";
      return 0;
    };

  private:

  };

} // algorithm

#endif
