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
 * @file lattice_shift_gpu.cc
 * GPU backend specialisation of Lattice_Shift_Tracker. Plain C++, no
 * CUDA language enabled (same split as leap_frog_gpu.cc/remove_com_
 * motion_gpu.cc); the real kernel lives in gpu/cuda/algorithm/
 * integration/lattice_shift_kernels.cu, reached through the host-safe
 * launch wrapper in lattice_shift_kernels.h.
 *
 * Only vacuum/rectangular boundary conditions are supported (matches
 * CUDA_Pairlist_Algorithm's own scope, cuda_pairlist_algorithm.cc) --
 * gpu::Periodicity<B>'s nearest_image() only implements those two
 * cases; init() hard-errors otherwise instead of silently producing
 * wrong box-wraps under triclinic/truncoct.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../util/template_split.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "../../gpu/cuda/algorithm/integration/lattice_shift_kernels.h"

#include "lattice_shift.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

#include "../../util/debug.h"

namespace {
  FPL9_TYPE matrix_to_fpl9(const math::Matrix & m) {
    FPL9_TYPE r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r(i, j) = static_cast<FPL_TYPE>(m(i, j));
    return r;
  }
}

template<>
int algorithm::Lattice_Shift_Tracker<util::gpuBackend>::
init(topology::Topology &topo,
     configuration::Configuration &conf,
     simulation::Simulation &sim,
     std::ostream &os,
     bool quiet) {

  const math::boundary_enum b = conf.boundary_type;
  if (b != math::vacuum && b != math::rectangular) {
    io::messages.add(
        "Lattice_Shift_Tracker<gpuBackend>: only vacuum and rectangular "
        "boundary conditions are supported so far (same scope as "
        "CUDA_Pairlist_Algorithm, cuda_pairlist_algorithm.cc); "
        "triclinic/truncoct is real future work, not yet done.",
        "Lattice_Shift_Tracker", io::message::error);
    return 1;
  }

  if (!sim.param().start.read_lattice_shifts) {
    conf.special().lattice_shifts = 0.0;
  }

  if (!quiet) {
    os << "LATTICESHIFTS (GPU)" << std::endl
       << "    keeping track of lattice shifts." << std::endl;

    if (sim.param().start.read_lattice_shifts)
      os << "    reading initial shifts from configuration.";
    else
      os << "    setting initial shifts to zero.";

    os << std::endl << "END" << std::endl;
  }

  return 0;
}

template<>
int algorithm::Lattice_Shift_Tracker<util::gpuBackend>::
apply(topology::Topology &topo,
      configuration::Configuration &conf,
      simulation::Simulation &sim) {
  DEBUG(6, "keeping track of lattice shifts (GPU)");
  this->m_timer.start(sim);
  SPLIT_BOUNDARY(_apply, topo, conf, sim);
  this->m_timer.stop();
  return 0;
}

template<>
template<math::boundary_enum b>
void algorithm::Lattice_Shift_Tracker<util::gpuBackend>::
_apply(topology::Topology &topo,
       configuration::Configuration &conf,
       simulation::Simulation &sim) {

  // Own stream: lets configuration_view() insert cudaStreamWaitEvent()
  // against whichever GPU-native algorithm produced POS last (typically
  // none yet this step -- this runs right after RemoveCOMMotion, near
  // the very top of the sequence) instead of relying on legacy-default-
  // stream ordering.
  static cudaStream_t stream = 0;
  if (stream == 0) cudaStreamCreate(&stream);

  gpu::Configuration::View view = sim.cuda().configuration_view(
      conf, gpu::MIRROR_POS | gpu::MIRROR_LATTICE_SHIFT, stream);
  const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

  // LATTICE_SHIFT has no "starts fresh" point anywhere else in the
  // codebase (deliberately excluded from MIRROR_ALL -- see mirror_
  // fields.h -- and this class's own gpu_mirror_touches() == 0 means
  // Algorithm_Sequence::run()'s generic invalidate never touches it
  // either): without this, mark_gpu_dirty() below would append one
  // more producer event to field_producer_events every single step,
  // forever, and every future configuration_view() request for this
  // field would cudaStreamWaitEvent() the entire ever-growing list --
  // an O(steps^2) cost that stayed invisible in short runs and only
  // dominated at 10000 steps (measured: this algorithm's own TIMING
  // line went 0.8s -> 17.5s over a full run before this fix). Only
  // clears the stale events (not the freshness bit, unlike invalidate_
  // gpu_mirror() -- this data stays GPU-resident, no CPU round trip):
  // configuration_view() above already resynced/waited on whatever was
  // pending, and mark_gpu_dirty() immediately below re-establishes
  // exactly one live event for this step's write.
  sim.cuda().clear_stale_producer_events(conf, gpu::MIRROR_LATTICE_SHIFT);

  // O(9) host math, exact CPU formula (math::Periodicity<b>::
  // put_chargegroups_into_box_saving_shifts(), math/periodicity.cc) --
  // cheap regardless of backend, not worth a kernel.
  const math::Box & my_box = conf.current().box;
  const math::Matrix L(my_box(0), my_box(1), my_box(2), true);
  const math::Matrix & cartesian_to_oblique = math::inverse(L);
  const FPL9_TYPE cartesian_to_oblique_fpl9 = matrix_to_fpl9(cartesian_to_oblique);

  gpu::launch_lattice_shift<b>(
      view.current().pos, view.lattice_shifts(),
      topo_view.chargegroup, topo_view.num_chargegroups,
      topo_view.num_solute_chargegroups,
      cartesian_to_oblique_fpl9, my_box, stream);

  // Vouch for what we just wrote: no CPU round trip. gpu_mirror_
  // touches() == 0 (lattice_shift.h) keeps the generic post-apply()
  // invalidation from immediately erasing this freshness.
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_POS | gpu::MIRROR_LATTICE_SHIFT, stream);
}

// explicit instantiations for linker
template class algorithm::Lattice_Shift_Tracker<util::gpuBackend>;
