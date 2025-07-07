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
 * @file gpu_settle.cc
 * contains the methods for the GPU_Settle class
 */


#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../math/periodicity.h"

#include "../../algorithm/constraints/shake.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"

#include "../../algorithm/constraints/gpu_settle.h"
// #ifdef HAVE_LIBCUDART
// #include "cukernel/cudaKernel.h"
// #endif

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

#define GPU_ID 0

int algorithm::GPU_Settle::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) {

  m_rank = 0;
  m_size = 1;


  if (!quiet) {
    os << "SETTLE\n"
       << "\tsolvent\n";
  }

  // check for 3 site water model like SPC, TIP3

  // we need a solvent
  if (topo.num_solvents() != 1) {
    io::messages.add("SETTLE does only work if 1 solvent.",
            "Settle", io::message::error);
    return 1;
  }

  // we need 3 atoms
  if (topo.solvent(0).num_atoms() != 3) {
    io::messages.add("SETTLE does only work with water like molecules (3 atoms).",
            "Settle", io::message::error);
    return 1;
  }

  // the masses of the second and third atom have to be the same
  if (topo.solvent(0).atom(1).mass != topo.solvent(0).atom(2).mass) {
    io::messages.add("SETTLE does only work with water like molecules (wrong masses).",
            "Settle", io::message::error);
    return 1;
  }

  // the molecule must be rigid: 3 distance constraints
  if (topo.solvent(0).distance_constraints().size() != 3) {
    io::messages.add("SETTLE does only work with water like molecules (3 distance constraints).",
            "Settle", io::message::error);
    return 1;
  }

  // the molecule must have two equal bond lengths (constraints 1 and 2)
  if (topo.bond_types_harm()[topo.solvent(0).distance_constraint(0).type].r0 !=
      topo.bond_types_harm()[topo.solvent(0).distance_constraint(1).type].r0) {
    io::messages.add("SETTLE does only work with water like molecules (distance constraints wrong).",
            "Settle", io::message::error);
    return 1;
  }

  // insert the solvent atoms to the constrained atoms
  for(unsigned int i = topo.num_solvent_atoms(); i < topo.num_atoms(); ++i) {
    constrained_atoms().insert(i);
  }

  // check whether we do an initial apply of the constraint algorithm
  if (sim.param().start.shake_pos || sim.param().start.shake_vel) {
    io::messages.add("initial settle-ing is not possible.",
            "Settle", io::message::error);
  }

  if (!quiet) {
    os << "END\n";
  }

  // gpu_stat = cukernel::cudaInitConstraints(sim.param().innerloop.number_gpus,
  //                                   sim.param().innerloop.gpu_device_number.at(GPU_ID),
  //                                   topo.num_solvent_atoms(0),
  //                                   topo.num_solvent_molecules(0));
  return 0;
}

/**
 * apply the SETTLE algorith
 */
int algorithm::GPU_Settle::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(7, "applying GPU_SETTLE");
  m_timer.start(sim);

  int error = 0;

  if (sim.param().system.nsm &&
      sim.param().constraint.solvent.algorithm == simulation::constr_gpu_settle) {

    DEBUG(8, "\twe need to settle SOLVENT");
    solvent(topo, conf, sim, error);

    if (error) {
      std::cout << "SETTLE: exiting with error condition "
              << "at step " << sim.steps() << std::endl;
      io::messages.add("SETTLE error", "Settle", io::message::error);
      conf.special().shake_failure_occurred = true;
      m_timer.stop();
      return 1;
    }
  }

  m_timer.stop();
  // return success!
  return 0;
}

/**
 * Let the CUDA device calculate the constraints
 */
void algorithm::GPU_Settle
::solvent(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim, int & error){

  // Backup the old positions
  math::Vec pos_bak[topo.num_solvent_atoms(0)];
  //pos_bak = (double *) malloc (sizeof(double) * topo.num_solvent_atoms(0));
  for (unsigned int i = 0; i < topo.num_solvent_atoms(0); ++i){
    pos_bak[i] = conf.current().pos(topo.num_solute_atoms() + i);
  }


  DEBUG(8, "Calculate constraints with GPU")
  int shake_fail_mol = -1;

  // Calculate the new poisitons on the GPU
  // error = cukernel::cudaConstraints(&conf.current().pos(topo.num_solute_atoms())(0),
  //                                         &conf.old().pos(topo.num_solute_atoms())(0),
  //                                         shake_fail_mol, gpu_stat);
  DEBUG(8, "GPU_Settle : Print Error, if one occured")
  if (error > 0)
    printError(topo, conf, sim, shake_fail_mol);

  // needed for constraint force, velocity correction and virial
  const double dt_i = 1.0 / sim.time_step_size();
  const double dt2_i = dt_i * dt_i;

  // masses
  const double mass_O = topo.solvent(0).atom(0).mass;
  assert(topo.solvent(0).atom(1).mass == topo.solvent(0).atom(2).mass);
  const double mass_H = topo.solvent(0).atom(1).mass;
  
  math::Vec cons_force [3];
  const bool do_velocity = !sim.param().stochastic.sd && !sim.param().minimise.ntem &&
          !sim.param().analyze.analyze;

  // Correct the velocities
  for (unsigned int i = 0; i < topo.num_solvent_atoms(); i += 3) {

    const math::Vec * const pos_old = &conf.old().pos(i);
    
    // calculate the displacement for velocity and virial calculation
    const math::Vec d_a = pos_bak[i + 0] - conf.current().pos(topo.num_solute_atoms() + i + 0);
    const math::Vec d_b = pos_bak[i + 1] - conf.current().pos(topo.num_solute_atoms() + i + 1);
    const math::Vec d_c = pos_bak[i + 2] - conf.current().pos(topo.num_solute_atoms() + i + 2);

    // in any case calculate the constraint force - it is a very interetsing
    // quantity
    // by finite difference
    cons_force[0] = d_a * mass_O * dt2_i;
    cons_force[1] = d_b * mass_H * dt2_i;
    cons_force[2] = d_c * mass_H * dt2_i;

    DEBUG(3, "constrained forces 0= " << math::v2s(cons_force[0]) << " 1= " << math::v2s(cons_force[1]) << " 2= " << math::v2s(cons_force[2]));

    if (do_velocity) {
      math::Vec * vel_new = &conf.current().vel(i);

      // let's reset the velocity
      // recalculate velocity from displacement and timestep
      vel_new[0] += d_a * dt_i;
      vel_new[1] += d_b * dt_i;
      vel_new[2] += d_c * dt_i;
    } // do velocity

    if (sim.param().pcouple.virial == math::atomic_virial) {
      // calculate the constraint virial
      for (int a = 0; a < 3; ++a) {
        for (int aa = 0; aa < 3; ++aa) {
          conf.old().virial_tensor(a, aa) +=
                  pos_old[0](a) * cons_force[0](aa) +
                  pos_old[1](a) * cons_force[1](aa) +
                  pos_old[2](a) * cons_force[2](aa);
        }
      }
    } // do virial
  } // for molecules
  return;
}

/**
 * a helper function to raise a settle error
 */
void algorithm::GPU_Settle
::printError(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        unsigned int mol) {
  unsigned int atom = topo.num_solute_atoms() + mol * topo.solvent(0).num_atoms();
  std::cout << "SETTLE ERROR\n"
          << "\tfailed to settle solvent molecule " << mol + 1 << ":\n"
          << "\tref OW:     " << math::v2s(conf.old().pos(atom)) << "\n"
          << "\tref H1:     " << math::v2s(conf.old().pos(atom + 1)) << "\n"
          << "\tref H2:     " << math::v2s(conf.old().pos(atom + 2)) << "\n"
          << "\tpos OW:     " << math::v2s(conf.current().pos(atom)) << "\n"
          << "\tpos H1:     " << math::v2s(conf.current().pos(atom + 1)) << "\n"
          << "\tpos H2:     " << math::v2s(conf.current().pos(atom + 2)) << "\n"
          << "\tvel OW:     " << math::v2s(conf.current().vel(atom)) << "\n"
          << "\tvel H1:     " << math::v2s(conf.current().vel(atom + 1)) << "\n"
          << "\tvel H2:     " << math::v2s(conf.current().vel(atom + 2)) << "\n"
          << "\told vel OW: " << math::v2s(conf.old().vel(atom)) << "\n"
          << "\told vel H1: " << math::v2s(conf.old().vel(atom + 1)) << "\n"
          << "\told vel H2: " << math::v2s(conf.old().vel(atom + 2)) << "\n"
          << "\tforce OW:   " << math::v2s(conf.old().force(atom)) << "\n"
          << "\tforce H1:   " << math::v2s(conf.old().force(atom + 1)) << "\n"
          << "\tforce H2:   " << math::v2s(conf.old().force(atom + 2)) << "\n";
}
