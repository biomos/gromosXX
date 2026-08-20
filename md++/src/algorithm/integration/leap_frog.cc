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
 * @file leap_frog.cc
 * CPU backend specialisations of Leap_Frog_Position and Leap_Frog_Velocity.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "leap_frog.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

// ─────────────────────────────────────────────────────────────────────────────
// Leap_Frog_Position<cpuBackend>
// x(t+dt) = x(t) + dt * v(t+dt/2)
// ─────────────────────────────────────────────────────────────────────────────

template <>
int algorithm::Leap_Frog_Position<util::cpuBackend>::apply(
    topology::Topology&           topo,
    configuration::Configuration& conf,
    simulation::Simulation&       sim)
{
    this->m_timer.start(sim);

    const int num_atoms = topo.num_atoms();
    const double dt     = sim.time_step_size();

#ifdef OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < num_atoms; ++i)
        conf.current().pos(i) = conf.old().pos(i) + conf.current().vel(i) * dt;

    if (sim.param().polarise.cos) {
#ifdef OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < num_atoms; ++i)
            conf.current().posV(i) = 2*conf.old().posV(i) - conf.current().posV(i);
    }

    this->m_timer.stop();
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Leap_Frog_Velocity<cpuBackend>
// v(t+dt/2) = v(t-dt/2) + dt/m * F(t)
// ─────────────────────────────────────────────────────────────────────────────

template <>
int algorithm::Leap_Frog_Velocity<util::cpuBackend>::apply(
    topology::Topology&           topo,
    configuration::Configuration& conf,
    simulation::Simulation&       sim)
{
    this->m_timer.start(sim);

    conf.exchange_state();
    sim.cuda().exchange_mirror_state(conf);
    conf.current().box = conf.old().box;

    const int    num_atoms = topo.num_atoms();
    const double dt        = sim.time_step_size();

#ifdef OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < num_atoms; ++i) {
        conf.current().vel(i) =
            conf.old().vel(i)
            + conf.old().force(i) * dt / topo.mass()(i);

        DEBUG(10, "atom " << i
              << "\n\tf=" << math::v2s(conf.old().force(i))
              << "\n\tmass=" << topo.mass()(i)
              << "\n\tvel=" << math::v2s(conf.old().vel(i)));
    }

    this->m_timer.stop();
    return 0;
}

// explicit instantiations for linker
template class algorithm::Leap_Frog_Position<util::cpuBackend>;
template class algorithm::Leap_Frog_Velocity<util::cpuBackend>;
