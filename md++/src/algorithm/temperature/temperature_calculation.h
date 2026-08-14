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
 * @file temperature/temperature_calculation.h
 * calculate the temperature.
 *
 * Backend-templated (PLAN.md §10 step 13 -- broadening algorithm
 * coverage). apply() is fully specialised per backend
 * (temperature_calculation_cpu.cc / temperature_calculation_gpu.cc);
 * init() has no backend-specific logic at all, so it stays a plain
 * inline template method here, same convention as
 * Leap_Frog_Velocity::init() (leap_frog.h).
 */

#ifndef INCLUDED_TEMPERATURE_CALCULATION_H
#define INCLUDED_TEMPERATURE_CALCULATION_H

#include "../../io/print_block.h"

namespace algorithm
{

  /**
   * @class Temperature_Calculation
   * temperature calculation.
   */
  template <typename Backend = util::cpuBackend>
  class Temperature_Calculation : public Algorithm, private AlgorithmB<Backend>
  {
  public:
    template <typename B>
    static constexpr bool is_supported_backend =
        std::is_same_v<B, util::cpuBackend> ||
        std::is_same_v<B, util::gpuBackend>;

    /**
     * Constructor.
     */
    Temperature_Calculation() : Algorithm("TemperatureCalculation") {}

    /**
     * Destructor.
     */
    virtual ~Temperature_Calculation() {}

    /**
     * apply the temperature calculation
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    /**
     * init -- no backend-specific logic, just prints the initial
     * multibath state after running apply() once.
     */
    virtual int init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false)
    {
      this->apply(topo, conf, sim);

      if (!quiet){
        io::print_MULTIBATH_COUPLING(os, sim.multibath());
        io::print_DEGREESOFFREEDOM(os, sim.multibath());
        io::print_MULTIBATH(os, sim.multibath(),
                            conf.old().energies,
                            "INITIAL TEMPERATURES");
      }
      return 0;
    }

    /**
     * gpuBackend never writes conf.current()/old() pos/vel directly --
     * only reads via the tracked accessor (sim.cuda().configuration_
     * view()), and only writes sim.multibath()/conf.old().energies,
     * which aren't part of the GPU-mirror freshness tracking at all.
     * Nothing to invalidate after this algorithm runs.
     */
    virtual unsigned gpu_mirror_touches() const override {
        if constexpr (std::is_same_v<Backend, util::gpuBackend>)
            return 0u;
        else
            return gpu::MIRROR_ALL;
    }

  private:

  };

} // algorithm

#endif
