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
 * @file nosehoover_thermostat.h
 * Nose-Hoover Thermostat.
 *
 * Backend-templated (PLAN.md §10 step 20), same shape as
 * Berendsen_Thermostat (berendsen_thermostat.h): calc_scaling()/
 * calc_chain_scaling() are O(num_baths) scalar math over
 * sim.multibath() -- never touch per-atom data, so they stay plain
 * inline template methods here. Only apply() (which does the
 * O(num_atoms) velocity scaling, via the inherited Thermostat::scale())
 * is specialised per backend (nosehoover_thermostat_cpu.cc /
 * nosehoover_thermostat_gpu.cc) -- and since Thermostat::scale()'s
 * uniform per-atom formula is exactly what Berendsen_Thermostat<
 * gpuBackend> already ported to a kernel (temperature_kernels.h's
 * launch_thermostat_scale_apply), NoseHoover_Thermostat<gpuBackend>
 * reuses that same kernel, not a new one.
 */

#ifndef INCLUDED_TEMPERATURE_NOSEHOOVER_H
#define INCLUDED_TEMPERATURE_NOSEHOOVER_H

#include "thermostat.h"

namespace algorithm
{

  /**
   * @class NoseHoover_Thermostat
   * the Nose-Hoover Thermostat.
   */
  template <typename Backend = util::cpuBackend>
  class NoseHoover_Thermostat : public Thermostat, private AlgorithmB<Backend>
  {
  public:
    template <typename B>
    static constexpr bool is_supported_backend =
        std::is_same_v<B, util::cpuBackend> ||
        std::is_same_v<B, util::gpuBackend>;

    /**
     * Constructor.
     */
    NoseHoover_Thermostat()
      : Thermostat("BerendsenThermostat") {}

    /**
     * Destructor.
     */
    virtual ~NoseHoover_Thermostat() {}

    /**
     * apply the temperature scaling
     * for baths with tau=-1 nothing is done.
     * the kinetic energy can not be calculated here, because
     * later on SHAKE might be applied.
     * @param topo the Topology
     * @param conf the Configuration
     * @param sim the Simulation
     */
    virtual int apply
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     );

    /**
     * initialise -- no backend-specific logic.
     */
    virtual int init
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     std::ostream & os = std::cout,
     bool quiet = false
     )
    {
      if (sim.param().multibath.algorithm > 0){

        if (!quiet){
          if (sim.param().multibath.algorithm == 1){
            std::cout << "\tNose-Hoover temperature coupling\n";
          }
          else{
            std::cout << "\tNose-Hoover-Chain temperature coupling: using "
                      << sim.param().multibath.algorithm << " instances\n";
          }
        }

        std::vector<simulation::bath_struct>::iterator
          b_it = sim.multibath().begin(),
          b_to = sim.multibath().end();

        for( ; b_it != b_to; ++b_it){
          b_it->zeta.resize(sim.param().multibath.algorithm, 0.0);
        }
      }
      return 0;
    }

    /**
     * calculate the scaling factors -- no backend-specific logic, pure
     * scalar math over sim.multibath(), never touches per-atom data.
     */
    void calc_scaling
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     )
    {
      std::vector<simulation::bath_struct>::iterator
        b_it = sim.multibath().begin(),
        b_to = sim.multibath().end();

      for(unsigned int num=0; b_it != b_to; ++b_it, ++num){
        if (b_it->tau != -1){

          double free_temp = 0.0;

          if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
            free_temp = 2 *
              (b_it->ekin - conf.special().flexible_constraint.flexible_ekin[num]) / (b_it->dof * math::k_Boltzmann);
          }
          else{
            free_temp = 2 *
              b_it->ekin / (b_it->dof * math::k_Boltzmann);
          }

          if (free_temp < math::epsilon) free_temp = b_it->temperature;

          b_it->zeta[0] += sim.time_step_size() / (b_it->tau * b_it->tau) * (free_temp / b_it->temperature - 1.0);
          b_it->scale = 1.0 - b_it->zeta[0] * sim.time_step_size();

        }
        else
          b_it->scale = 1;
      }
    }

    /**
     * calculate the scaling factors for Nose-Hoover Chains -- no
     * backend-specific logic, pure scalar math over sim.multibath().
     */
    void calc_chain_scaling
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     )
    {
      const double dt = sim.time_step_size();

      std::vector<simulation::bath_struct>::iterator
        b_it = sim.multibath().begin(),
        b_to = sim.multibath().end();

      for(unsigned int num=0; b_it != b_to; ++b_it, ++num){
        if (b_it->tau != -1){

          double free_temp = 0.0;

          if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
            free_temp = 2 *
              (b_it->ekin - conf.special().flexible_constraint.flexible_ekin[num]) /
              (b_it->dof * math::k_Boltzmann);
          }
          else{
            free_temp = 2 *
              b_it->ekin / (b_it->dof * math::k_Boltzmann);
          }

          if (free_temp < math::epsilon) free_temp = b_it->temperature;

          const int nhc = sim.param().multibath.algorithm;

          std::vector<double> tau(nhc, b_it->tau * b_it->tau / b_it->dof);
          tau[0] = b_it->tau * b_it->tau;

          assert(nhc > 1);

          b_it->zeta[nhc-1] += (tau[nhc-2] * b_it->zeta[nhc-2] * b_it->zeta[nhc-2]
                             - 1.0 / b_it->dof) / tau[nhc-1] * dt;

          for (int i = nhc - 2; i >= 1; i--){

            b_it->zeta[i] += ((tau[i-1] * b_it->zeta[i-1] * b_it->zeta[i-1] - 1.0 / b_it->dof)
                             / tau[i]  - b_it->zeta[i] * b_it->zeta[i+1]) * dt;
          }

          b_it->zeta[0] += ((free_temp / b_it->temperature - 1.0) / tau[0]
                             - b_it->zeta[1] * b_it->zeta[0] ) * dt;

          b_it->scale = 1.0 - b_it->zeta[0] * dt;

        }
        else
          b_it->scale = 1;
      }
    }

    /**
     * gpuBackend leaves the scaled velocity resident on the GPU mirror
     * -- same reasoning as Berendsen_Thermostat<gpuBackend>.
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
