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
 * @file cuda_molecular_virial_interaction.cc
 * GPU-native molecular-virial correction. See
 * cuda_molecular_virial_interaction.h. Plain C++, no kernel syntax; the
 * real __global__ kernels live in
 * gpu/cuda/interaction/special/molecular_virial_kernels.cu, compiled
 * into grocuda.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "gpu/cuda/interaction/special/molecular_virial_kernels.h"
#include "gpu/cuda/memory/virial_accumulate_kernels.h"

#include "cuda_molecular_virial_interaction.h"

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

int interaction::CUDA_Molecular_Virial_Interaction::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  const std::vector<unsigned int> & pg = topo.pressure_groups();
  m_num_groups = static_cast<unsigned>(pg.size()) - 1;
  m_num_atoms = static_cast<unsigned>(topo.num_atoms());

  m_group_offsets.resize(pg.size());
  for (unsigned i = 0; i < pg.size(); ++i) {
    m_group_offsets[i] = pg[i];
  }

  m_group_id.resize(m_num_atoms);
  for (unsigned g = 0; g < m_num_groups; ++g) {
    for (unsigned a = pg[g]; a < pg[g + 1]; ++a) {
      m_group_id[a] = g;
    }
  }

  m_com_pos.resize(m_num_groups);
  m_virial.resize(9);

  if (m_stream == 0) cudaStreamCreate(&m_stream);

  m_initialized = true;

  if (!quiet)
    os << "CUDA MOLECULAR VIRIAL INTERACTION\n"
       << "\tpressure groups: " << m_num_groups << "\n"
       << "END\n";
  return 0;
}

interaction::CUDA_Molecular_Virial_Interaction::~CUDA_Molecular_Virial_Interaction() {
  if (m_stream) cudaStreamDestroy(m_stream);
}

int interaction::CUDA_Molecular_Virial_Interaction::calculate_interactions(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  cudaMemsetAsync(m_virial.data(), 0, sizeof(double) * 9, m_stream);

  // Needs both POS (for the per-group unwrap and the per-atom nearest-
  // image to the group COM) and FORCE (the accumulated atomic force
  // bonded/nonbonded terms already wrote into the mirror this step) --
  // both read-only here, never written back to the mirror by this
  // class.
  gpu::Configuration::View view =
      sim.cuda().configuration_view(conf, gpu::MIRROR_POS | gpu::MIRROR_FORCE, m_stream);
  const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

  gpu::launch_group_com(
      view.current().pos, topo_view.mass, m_group_offsets.data(), m_num_groups,
      conf.boundary_type, conf.current().box, m_com_pos.data(), m_stream);

  // Same stream as launch_group_com() above -- CUDA's single-stream
  // in-order execution guarantee is what makes this pass safe to read
  // m_com_pos without an explicit host-side sync between the two.
  gpu::launch_molecular_virial_correction(
      view.current().pos, view.current().force.data(), m_group_id.data(), m_com_pos.data(),
      m_num_atoms, conf.boundary_type, conf.current().box, m_virial.data(), m_stream);

  // Publish into the shared mirror's virial_tensor via atomicAdd
  // (GPU-resident, no CPU round trip) -- m_virial already holds the
  // negated correction (see molecular_virial_kernels.h's doc comment),
  // so this atomicAdd is exactly virial_tensor -= corrP.
  gpu::launch_accumulate_virial9(
      reinterpret_cast<FPH_TYPE*>(view.current().virial_tensor), m_virial.data(), m_stream);
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VIRIAL, m_stream);

  // Nothing to read back to the host -- no cudaStreamSynchronize() at
  // all, unlike the bonded-term classes (which need one for their
  // private per-energy-group buffer). This class has no energy output.
  m_timer.stop();
  return 0;
}
