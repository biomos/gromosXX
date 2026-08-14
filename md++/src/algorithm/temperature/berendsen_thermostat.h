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
 * @file berendsen_thermostat.h
 * berendsen thermostat.
 *
 * Backend-templated (PLAN.md §10 step 13). calc_scaling() is O(num_baths)
 * scalar math over sim.multibath() -- never touches per-atom data, so it
 * stays a plain inline template method here (no backend-specific logic
 * at all), same convention as Temperature_Calculation::init(). Only
 * apply() (which does the O(num_atoms) velocity scaling) is fully
 * specialised per backend (berendsen_thermostat_cpu.cc /
 * berendsen_thermostat_gpu.cc).
 */

#ifndef INCLUDED_TEMPERATURE_BERENDSEN_H
#define INCLUDED_TEMPERATURE_BERENDSEN_H

#include "thermostat.h"

namespace algorithm
{

  /**
   * @class Berendsen_Thermostat
   * the Berendsen thermostat.
   */
  template <typename Backend = util::cpuBackend>
  class Berendsen_Thermostat : public Thermostat, private AlgorithmB<Backend>
  {
  public:
    template <typename B>
    static constexpr bool is_supported_backend =
        std::is_same_v<B, util::cpuBackend> ||
        std::is_same_v<B, util::gpuBackend>;

    /**
     * Constructor.
     */
    Berendsen_Thermostat() : Thermostat("BerendsenThermostat") {}

    /**
     * Destructor.
     */
    virtual ~Berendsen_Thermostat() {}

    /**
     * initialise -- no backend-specific logic.
     */
    virtual int init(topology::Topology & topo,
                     configuration::Configuration & conf,
                     simulation::Simulation & sim,
                     std::ostream & os = std::cout,
                     bool quiet = false)
    {
      if (!quiet){
        os << "\tWeak-Coupling temperature coupling\n";
      }
      return 0;
    }

    /**
     * apply the temperature scaling
     * for baths with tau=-1 nothing is done.
     * the kinetic energy can not be calculated here, because
     * later on SHAKE might be applied.
     * @param topo the Topology
     * @param conf the Configuration
     * @param sim the Simulation
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    /**
     * calculate the scaling factors -- no backend-specific logic, pure
     * scalar math over sim.multibath(), never touches per-atom data.
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param immediate if true rescales the velocities to immediately satisfy
     * the given reference temperature (strong coupling).
     */
    void calc_scaling(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      bool immediate = false)
    {
      std::vector<simulation::bath_struct>::iterator
        b_it = sim.multibath().begin(),
        b_to = sim.multibath().end();

      for(unsigned int num=0; b_it != b_to; ++b_it, ++num){
        if (b_it->tau != -1 || immediate){

          double free_temp = 0.0;

          if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
            free_temp = 2 *
              (b_it->ekin - conf.special().flexible_constraint.flexible_ekin[num])
              / (b_it->dof * math::k_Boltzmann);
          }
          else{
            free_temp = 2 *
              b_it->ekin / (b_it->dof * math::k_Boltzmann);
          }

          if (free_temp < math::epsilon) free_temp = b_it->temperature;

          if (free_temp < math::epsilon)
            b_it->scale = 1;
          else if (immediate)
            b_it->scale = sqrt(b_it->temperature / free_temp);
          else
            b_it->scale = sqrt(1.0 + sim.time_step_size() / b_it->tau *
                               (b_it->temperature / free_temp - 1));
        }
        else
          b_it->scale = 1;
      }
    }

    /**
     * gpuBackend leaves the scaled velocity resident on the GPU mirror
     * (mark_gpu_dirty(MIRROR_VEL), no sync-back) -- same reasoning as
     * Leap_Frog_Velocity<gpuBackend>: this sits between Leap_Frog_
     * Velocity and Leap_Frog_Position in create_md_sequence.cc, so
     * returning 0 here keeps Algorithm_Sequence::run()'s default
     * post-apply() invalidation from immediately erasing that freshness,
     * letting Leap_Frog_Position read the scaled velocity with zero
     * extra round trip.
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
