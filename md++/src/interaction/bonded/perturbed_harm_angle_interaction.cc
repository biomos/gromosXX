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
 * @file perturbed_harm_angle_interaction.cc
 * methods of Perturbed_Harm_Angle_Interaction.
 */

#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"
#include "interaction/interaction.h"
#include "math/periodicity.h"
#include "perturbed_harm_angle_interaction.h"
#include "util/template_split.h"

template<math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_harm_angle_interactions
(topology::Topology &topo,
 configuration::Configuration &conf,
 simulation::Simulation &sim,
 int &error)
{
  std::vector<interaction::angle_type_struct> const &param = topo.angle_types_harm();
  std::vector<topology::perturbed_three_body_term_struct>::const_iterator a_it =
    topo.perturbed_solute().angles().begin(),
    a_to = topo.perturbed_solute().angles().end();

  math::VArray &pos = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Periodicity<B> periodicity(conf.current().box);

  for( ; a_it != a_to; ++a_it){
    const unsigned int group = topo.atom_energy_group()[a_it->i];
    const double lambda = topo.individual_lambda(simulation::angle_lambda)[group][group];
    const double lambda_derivative =
      topo.individual_lambda_derivative(simulation::angle_lambda)[group][group];

    assert(unsigned(a_it->A_type) < param.size());
    assert(unsigned(a_it->B_type) < param.size());

    const double K = (1 - lambda) * param[a_it->A_type].K + lambda * param[a_it->B_type].K;
    // cos0 contains an angle in radians in the harmonic parameter table.
    const double theta0 = (1 - lambda) * param[a_it->A_type].cos0 +
      lambda * param[a_it->B_type].cos0;
    const double K_diff = param[a_it->B_type].K - param[a_it->A_type].K;
    const double theta_diff = param[a_it->B_type].cos0 - param[a_it->A_type].cos0;

    math::Vec rij, rkj;
    periodicity.nearest_image(pos(a_it->i), pos(a_it->j), rij);
    periodicity.nearest_image(pos(a_it->k), pos(a_it->j), rkj);
    const double dij = sqrt(abs2(rij));
    const double dkj = sqrt(abs2(rkj));
    if (!(dij > 0.0 && dkj > 0.0) || !std::isfinite(dij) || !std::isfinite(dkj)){
      io::messages.add("zero or nonfinite bond length in harmonic angle",
                      "Perturbed_Harm_Angle_Interaction", io::message::critical);
      error = 1;
      continue;
    }

    double cost = dot(rij, rkj) / (dij * dkj);
    if (!std::isfinite(cost)){
      io::messages.add("nonfinite harmonic angle geometry",
                      "Perturbed_Harm_Angle_Interaction", io::message::critical);
      error = 1;
      continue;
    }
    cost = std::max(-1.0, std::min(1.0, cost));
    const double theta = acos(cost);
    const double sint = sin(theta);
    const double diff = theta - theta0;
    double energy = 0.5 * K * diff * diff;
    double e_lambda = lambda_derivative * (0.5 * K_diff * diff * diff - K * diff * theta_diff);
    double fac;
    if (sint < math::epsilon){
      if (cost > 0.0 || fabs(theta0 - math::Pi) > math::epsilon){
        io::messages.add("singular harmonic angle away from 180 degree equilibrium",
                        "Perturbed_Harm_Angle_Interaction", io::message::critical);
        error = 1;
        continue;
      }
      // Same theta -> pi limit as Harm_Angle_Interaction.
      fac = -K;
    }
    else{
      fac = K * diff / sint;
    }

    const math::Vec fi = fac / dij * (rkj / dkj - rij / dij * cost);
    const math::Vec fk = fac / dkj * (rij / dij - rkj / dkj * cost);
    const math::Vec fj = -fi - fk;
    force(a_it->i) += fi;
    force(a_it->j) += fj;
    force(a_it->k) += fk;

    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
        conf.current().virial_tensor(a, b) += rij(a) * fi(b) + rkj(a) * fk(b);

    conf.current().energies.angle_energy[group] += energy;
    conf.current().perturbed_energy_derivatives.angle_energy[group] += e_lambda;

    if (sim.param().gamd.gamd){
      const unsigned int gamd_group = topo.gamd_accel_group(a_it->i);
      std::vector<unsigned int> key = {gamd_group, gamd_group};
      const unsigned int igroup = topo.gamd_interaction_group(key);
      conf.special().gamd.total_force[igroup](a_it->i) += fi;
      conf.special().gamd.total_force[igroup](a_it->j) += fj;
      conf.special().gamd.total_force[igroup](a_it->k) += fk;
      conf.current().energies.gamd_potential_total[igroup] += energy;
      for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
          conf.special().gamd.virial_tensor[igroup](a, b) += rij(a) * fi(b) + rkj(a) * fk(b);
    }

    if (sim.param().precalclam.nr_lambdas &&
        sim.param().write.free_energy > 0 &&
        sim.steps() % sim.param().write.free_energy == 0){
      const double lambda_step = sim.param().precalclam.nr_lambdas > 1 ?
        (sim.param().precalclam.max_lam - sim.param().precalclam.min_lam) /
        (sim.param().precalclam.nr_lambdas - 1) : 0.0;
      for(unsigned int l=0; l<sim.param().precalclam.nr_lambdas; ++l){
        const double lam = sim.param().precalclam.min_lam + l * lambda_step;
        const double Klam = (1 - lam) * param[a_it->A_type].K + lam * param[a_it->B_type].K;
        const double theta0lam = (1 - lam) * param[a_it->A_type].cos0 + lam * param[a_it->B_type].cos0;
        const double difflam = theta - theta0lam;
        conf.current().energies.AB_angle[l] += 0.5 * Klam * difflam * difflam;
        conf.current().perturbed_energy_derivatives.AB_angle[l] +=
          0.5 * K_diff * difflam * difflam - Klam * difflam * theta_diff;
      }
    }
  }
  return 0;
}

int interaction::Perturbed_Harm_Angle_Interaction
::calculate_interactions(topology::Topology &topo,
                         configuration::Configuration &conf,
                         simulation::Simulation &sim)
{
  int error = 0;
  m_timer.start(sim);
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_harm_angle_interactions, topo, conf, sim, error);
  m_timer.stop();
  return error;
}
