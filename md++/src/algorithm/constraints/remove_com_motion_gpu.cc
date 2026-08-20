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
 * @file remove_com_motion_gpu.cc
 * remove com motion - gpu variant. Plain C++, no kernel syntax
 * (groalgorithm has no CUDA language enabled -- same constraint
 * leap_frog_gpu.cc already works within); the real __global__ kernels
 * live in gpu/cuda/algorithm/constraints/remove_com_motion_kernels.cu,
 * compiled into grocuda, reached only through the host-safe launch
 * wrappers declared in remove_com_motion_kernels.h.
 *
 * Mirrors remove_com_motion_cpu.cc's math exactly (same formulas, same
 * order of operations) -- see that file for the derivation. This one
 * does the reductions (sum(m*v), sum(m*pos - 0.5*m*v*dt), angular
 * momentum, inertia tensor) on the GPU instead of a CPU loop; the tiny
 * follow-up scalar work (dividing by total mass, inverting the 3x3
 * inertia tensor) stays on the host -- not worth a kernel for 9
 * numbers.
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "../../gpu/cuda/algorithm/constraints/remove_com_motion_kernels.h"

#include "remove_com_motion.h"

#include "../../io/print_block.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

template<>
int algorithm::Remove_COM_Motion<util::gpuBackend>::init
(
 topology::Topology &topo,
 configuration::Configuration &conf,
 simulation::Simulation &sim,
 std::ostream &os,
 bool quiet
)
{
  if (quiet) return 0;

  os << "CENTRE OF MASS MOTION (GPU)\n";

  if (sim.param().centreofmass.skip_step){
    if (sim.param().centreofmass.skip_step > 1)
      os << "\tremoving centre of mass motion every "
	 << sim.param().centreofmass.skip_step
	 << " steps\n";
    else
      os << "\tremoving centre of mass motion every step\n";

    if (sim.param().centreofmass.remove_rot){
      os << "\tremoving centre of mass rotation" << std::endl;
    }
    if (sim.param().centreofmass.remove_trans){
      os << "\tremoving centre of mass translation" << std::endl;
    }
    os << "\n";
  }
  if (sim.param().print.centreofmass > 1){
    os << "\tprinting centre of mass motion every "
       << sim.param().print.centreofmass
       << " steps\n";
  }
  if (sim.param().print.centreofmass == 1){
    os << "\tprinting centre of mass motion every step\n";
  }

  if (sim.param().start.remove_com_translation)
    os << "\n\tremoving initial centre of mass translation\n";
  if (sim.param().start.remove_com_rotation)
    os << "\n\tremoving initial centre of mass rotation\n";

  os << "END\n";

  return 0;
}

template<>
double algorithm::Remove_COM_Motion<util::gpuBackend>
::remove_com_translation
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 bool remove_trans
 )
{
  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());

  // Own stream, function-local static (not a class member) since
  // Remove_COM_Motion's header compiles unconditionally in CPU-only
  // builds too and can't hold a cudaStream_t directly. This is a
  // stream-scoped cudaStreamSynchronize() below, not a device-wide
  // cudaDeviceSynchronize(), so it doesn't stall whatever else (bonded
  // terms, NonBonded, constraints) is running concurrently on its own
  // stream -- still a genuine CPU-blocking wait though: the reduction's
  // result (com_v_x/y/z) is a host-computed launch parameter for
  // launch_com_translation_apply() right below, an inherent sequential
  // dependency this multi-pass-reduction design can't route around via
  // the event mechanism (that only helps GPU-GPU ordering).
  static cudaStream_t stream = 0;
  if (stream == 0) cudaStreamCreate(&stream);

  gpu::Configuration::View conf_view =
      sim.cuda().configuration_view(conf, gpu::MIRROR_VEL, stream);
  const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

  // Persistent, GPU-only scratch for the reduction -- fixed size (4
  // doubles), allocated once on first call and reused thereafter.
  static gpu::cuvector<double> sums;
  if (sums.size() < 4) sums.resize(4);

  gpu::launch_com_translation_reduce(conf_view.current().vel, topo_view.mass, num_atoms, sums.data(), stream);
  cudaStreamSynchronize(stream);

  const double com_mass = sums[3];
  const double com_v_x = sums[0] / com_mass;
  const double com_v_y = sums[1] / com_mass;
  const double com_v_z = sums[2] / com_mass;

  const double ekin_trans = 0.5 * com_mass *
      (com_v_x * com_v_x + com_v_y * com_v_y + com_v_z * com_v_z);

  if (remove_trans) {
    gpu::launch_com_translation_apply(conf_view.current().vel, com_v_x, com_v_y, com_v_z, num_atoms, stream);
    sim.cuda().sync_configuration_from_device(conf);
  }

  return ekin_trans;
}

template<>
double algorithm::Remove_COM_Motion<util::gpuBackend>
::remove_com_rotation
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 bool remove_rot
 )
{
  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  const double dt = sim.time_step_size();

  // See remove_com_translation()'s comment: own stream, but the two
  // reduction passes are an inherent host-in-the-loop sequential
  // dependency (pass2's launch parameters are pass1's host-side
  // result), not something the event mechanism can route around.
  static cudaStream_t stream = 0;
  if (stream == 0) cudaStreamCreate(&stream);

  gpu::Configuration::View conf_view =
      sim.cuda().configuration_view(conf, gpu::MIRROR_POS | gpu::MIRROR_VEL, stream);
  const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

  static gpu::cuvector<double> sums1;
  static gpu::cuvector<double> sums2;
  if (sums1.size() < 7)  sums1.resize(7);
  if (sums2.size() < 12) sums2.resize(12);

  gpu::launch_com_rotation_reduce_pass1(
      conf_view.current().pos, conf_view.current().vel, topo_view.mass, dt, num_atoms, sums1.data(), stream);
  cudaStreamSynchronize(stream);

  const double com_mass = sums1[6];
  const double com_v_x = sums1[0] / com_mass, com_v_y = sums1[1] / com_mass, com_v_z = sums1[2] / com_mass;
  const double com_r_x = sums1[3] / com_mass, com_r_y = sums1[4] / com_mass, com_r_z = sums1[5] / com_mass;

  DEBUG(7, "com_r " << com_r_x << " " << com_r_y << " " << com_r_z);
  DEBUG(7, "com_v " << com_v_x << " " << com_v_y << " " << com_v_z);

  gpu::launch_com_rotation_reduce_pass2(
      conf_view.current().pos, conf_view.current().vel, topo_view.mass, dt,
      com_v_x, com_v_y, com_v_z, com_r_x, com_r_y, com_r_z, num_atoms, sums2.data(), stream);
  cudaStreamSynchronize(stream);

  const double Lx = sums2[0], Ly = sums2[1], Lz = sums2[2];
  math::Matrix com_I;
  com_I(0,0) = sums2[3];  com_I(0,1) = sums2[4];  com_I(0,2) = sums2[5];
  com_I(1,0) = sums2[6];  com_I(1,1) = sums2[7];  com_I(1,2) = sums2[8];
  com_I(2,0) = sums2[9];  com_I(2,1) = sums2[10]; com_I(2,2) = sums2[11];

  DEBUG(7, "Angular momentum " << Lx << " " << Ly << " " << Lz);

  // Exact CPU formula (remove_com_motion_cpu.cc) -- invert the inertia
  // tensor, trivial 3x3 work, not worth a kernel for 9 numbers.
  math::Matrix com_II;
  const double denom = -com_I(2,0)*com_I(2,0)*com_I(1,1)
    + 2 * com_I(0,1) * com_I(0,2) * com_I(1,2)
    - com_I(0, 0) * com_I(1,2) * com_I(1,2)
    - com_I(0,1) * com_I(0,1) * com_I(2,2)
    + com_I(0,0) * com_I(1,1) * com_I(2,2);

  com_II(0,0) = (-com_I(1,2)*com_I(1,2) + com_I(1,1) * com_I(2,2));
  com_II(1,0) = com_II(0,1) = (com_I(0,2) * com_I(1,2)
    - com_I(0,1) * com_I(2,2));
  com_II(0,2) = com_II(2,0) = (-com_I(0,2)*com_I(1,1)
    + com_I(0,1)*com_I(1,2));

  com_II(1,1) = (-com_I(0,2)*com_I(0,2) + com_I(0,0) * com_I(2,2));
  com_II(1,2) = com_II(2,1) = (com_I(0,1)*com_I(0,2)
    - com_I(0,0) * com_I(1,2));

  com_II(2,2) = (-com_I(0,1)*com_I(0,1) + com_I(0,0)*com_I(1,1));

  DEBUG(7, "inertia tensor:\n"<< math::m2s(com_I));
  DEBUG(7, "determinant : " << denom);
  DEBUG(7, "inverted tens :\n" << math::m2s(com_II));

  math::Vec com_L(Lx, Ly, Lz);
  math::Vec com_O(0.0);
  if (denom >= math::epsilon)
    com_O = math::product(com_II, com_L) / denom;

  DEBUG(7, " angular velocity " << math::v2s(com_O));

  const double ekin_rot = 0.5 * dot(com_O, com_L);
  DEBUG(7, " com_Ekin_rot " << ekin_rot);

  if (remove_rot) {
    gpu::launch_com_rotation_apply(
        conf_view.current().pos, conf_view.current().vel, dt,
        com_r_x, com_r_y, com_r_z, com_O(0), com_O(1), com_O(2), num_atoms, stream);
    sim.cuda().sync_configuration_from_device(conf);
  }

  return ekin_rot;
}

/**
 * apply the COM removal.
 */
template<>
int algorithm::Remove_COM_Motion<util::gpuBackend>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  bool remove_rot = false;
  bool remove_trans = false;
  bool print_it = false;

  m_timer.start(sim);

  // check if nothing to do
  if (sim.steps() == 0){
    remove_rot = sim.param().start.remove_com_rotation;
    remove_trans = sim.param().start.remove_com_translation;
    if (sim.param().print.centreofmass) print_it = true;
  }
  else{
    if (sim.param().centreofmass.skip_step &&
	(sim.steps() % sim.param().centreofmass.skip_step) == 0)
      remove_rot = remove_trans = true;
    if (sim.param().print.centreofmass &&
	(sim.steps() % sim.param().print.centreofmass) == 0)
      print_it = true;
  }

  DEBUG(9, "centre of mass: print " << print_it << " remove " <<
           (remove_trans || remove_rot) );
  if (!print_it && !remove_trans && !remove_rot){
    m_timer.stop();
    return 0;
  }

  if (sim.steps() != 0){
    remove_rot = remove_rot && sim.param().centreofmass.remove_rot;
    remove_trans = remove_trans && sim.param().centreofmass.remove_trans;
  }

  DEBUG(9, "centre of mass: trans " << remove_trans << " rot " << remove_rot);

  double ekin_trans = 0.0, ekin_rot = 0.0;

  if (print_it || remove_trans){
    ekin_trans = remove_com_translation(topo, conf, sim, remove_trans);
  }
  if (print_it || remove_rot){
    ekin_rot = remove_com_rotation(topo, conf, sim, remove_rot);
  }

  if (print_it){
    io::print_CENTREOFMASS(os, ekin_trans, ekin_rot);
  }

  m_timer.stop();

  return 0;
}

template<>
double algorithm::Remove_COM_Motion<util::gpuBackend>
::add_com_rotation
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 math::Vec com_L
 )
{
  // Never called anywhere in this codebase (confirmed by grep) --
  // no GPU implementation, hard-error rather than silently doing
  // nothing if that ever changes.
  io::messages.add(
      "Remove_COM_Motion<gpuBackend>::add_com_rotation is not implemented "
      "(unused in this codebase today -- see remove_com_motion_gpu.cc).",
      "Remove_COM_Motion", io::message::error);
  return 0.0;
}

// explicit instantiations for linker
template class algorithm::Remove_COM_Motion<util::gpuBackend>;
