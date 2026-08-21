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
 * @file berendsen_barostat_gpu.cc
 * GPU backend specialisation of Berendsen_Barostat. Plain C++ -- no CUDA
 * language enabled here (same split as leap_frog_gpu.cc/remove_com_
 * motion_gpu.cc); the real __global__ kernel lives in gpu/cuda/
 * algorithm/pressure/berendsen_barostat_kernels.cu, compiled into
 * grocuda, reached only through the host-safe launch wrapper declared
 * in berendsen_barostat_kernels.h.
 *
 * mu (the 3x3 scaling matrix) and the box update are computed on the
 * host exactly like Berendsen_Barostat<cpuBackend>::apply()
 * (berendsen_barostat_cpu.cc) -- both are O(9), trivially cheap
 * regardless of backend, and box is plain CPU-resident data. Only
 * position scaling (O(num_atoms)) is GPU-resident: one generic
 * matrix-vector kernel handles isotropic/anisotropic/semi-anisotropic
 * alike, since all three reduce to "multiply every position by a fixed
 * (here: diagonal) 3x3 matrix." Positions stay on the GPU mirror across
 * this call -- no eager publish; see gpu_mirror_touches() in
 * berendsen_barostat.h.
 *
 * full_anisotropic is deliberately NOT ported: reading berendsen_
 * barostat_cpu.cc closely, its full_anisotropic case block has no
 * `break;` before the following `case math::pcouple_semi_anisotropic:`
 * -- a pre-existing fallthrough bug in the CPU reference that makes
 * full_anisotropic mode also re-run (and re-apply) the semi-anisotropic
 * scaling on top of its own result. Silently reproducing that on GPU
 * would enshrine it as "the answer"; silently fixing it here would
 * change CPU-observable behaviour as an undiscussed side effect of a
 * performance patch. Neither is this task's call to make, so gpuBackend
 * refuses this one mode with a clear error instead (see apply() below)
 * -- not exercised by any benchmark/test in this branch.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../configuration/state_properties.h"

#include "../../util/error.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "../../gpu/cuda/algorithm/pressure/berendsen_barostat_kernels.h"

#include "berendsen_barostat.h"
#include "../../math/boundary_checks.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE pressure

#include "../../util/debug.h"
#include "../../math/transformation.h"

namespace {
  FPL9_TYPE diag_to_fpl9(const math::Vec & d) {
    FPL9_TYPE r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r(i, j) = static_cast<FPL_TYPE>(i == j ? d(i) : 0.0);
    return r;
  }
}

template<>
int algorithm::Berendsen_Barostat<util::gpuBackend>
::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  m_timer.start(sim);

  DEBUG(8, "Berendsen Barostat (GPU) == apply");

  math::Matrix & pressure = conf.old().pressure_tensor;
  math::Box & box = conf.current().box;

  bool scale_ref = sim.param().posrest.posrest != simulation::posrest_off &&
          sim.param().posrest.scale_reference_positions;

  // mu(0), mu(1), mu(2): the diagonal of the per-atom scaling matrix
  // (isotropic: all three equal; anisotropic/semi-anisotropic: per-axis).
  math::Vec mu;

  switch (sim.param().pcouple.scale) {
    case math::pcouple_isotropic:
    {
      double total_pressure = (pressure(0, 0) + pressure(1, 1) + pressure(2, 2)) / 3.0;
      const double m = pow(1.0 - sim.param().pcouple.compressibility
              * sim.time_step_size() / sim.param().pcouple.tau
              * (sim.param().pcouple.pres0(0, 0) - total_pressure),
              1.0 / 3.0);
      mu = math::Vec(m, m, m);
      box *= m;
      break;
    }
    case math::pcouple_anisotropic:
    {
      for (int i = 0; i < 3; ++i) {
        mu(i) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * (sim.param().pcouple.pres0(i, i) - pressure(i, i)),
                1.0 / 3.0);
      }
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          box(i)(j) *= mu(j);
      break;
    }
    case math::pcouple_semi_anisotropic:
    {
      if (sim.param().pcouple.x_semi < 1)
        mu(0) = 1;
      else {
        mu(0) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * (sim.param().pcouple.pres0(0, 0) - pressure(0, 0)),
                1.0 / 3.0);
      }

      if (sim.param().pcouple.y_semi < 1)
        mu(1) = 1;
      if (sim.param().pcouple.y_semi > 0 && (sim.param().pcouple.y_semi == sim.param().pcouple.x_semi)) {
        mu(0) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * ((sim.param().pcouple.pres0(0, 0) + sim.param().pcouple.pres0(1, 1) -
                pressure(0, 0) - pressure(1, 1)) / 2),
                1.0 / 3.0);
        mu(1) = mu(0);
      }
      if (sim.param().pcouple.y_semi > 0 && (sim.param().pcouple.y_semi != sim.param().pcouple.x_semi)) {
        mu(1) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * (sim.param().pcouple.pres0(1, 1) - pressure(1, 1)),
                1.0 / 3.0);
      }

      if (sim.param().pcouple.z_semi < 1)
        mu(2) = 1;
      if (sim.param().pcouple.z_semi > 0) {
        if (sim.param().pcouple.z_semi == sim.param().pcouple.x_semi) {
          mu(0) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * ((sim.param().pcouple.pres0(0, 0) + sim.param().pcouple.pres0(2, 2) -
                pressure(0, 0) - pressure(2, 2)) / 2),
                1.0 / 3.0);
          mu(2) = mu(0);
        }
        if (sim.param().pcouple.z_semi == sim.param().pcouple.y_semi) {
          mu(1) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * ((sim.param().pcouple.pres0(1, 1) + sim.param().pcouple.pres0(2, 2) -
                pressure(1, 1) - pressure(2, 2)) / 2),
                1.0 / 3.0);
          mu(2) = mu(1);
        }
        if ((sim.param().pcouple.z_semi != sim.param().pcouple.x_semi) &&
            (sim.param().pcouple.z_semi != sim.param().pcouple.y_semi)) {
          mu(2) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * (sim.param().pcouple.pres0(2, 2) - pressure(2, 2)),
                1.0 / 3.0);
        }
      }

      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          box(i)(j) *= mu(j);
      break;
    }
    case math::pcouple_full_anisotropic:
    {
      io::messages.add(
          "Berendsen_Barostat<gpuBackend>: full-anisotropic pressure "
          "scaling is not supported on the GPU integration path.",
          "Berendsen_Barostat", io::message::error);
      m_timer.stop();
      return 1;
    }
    default:
      m_timer.stop();
      return 0;
  }

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());

  // Reference positions (position-restraint scaling) are never part of
  // the GPU mirror (see cuda_position_restraint_interaction.h's doc
  // comment: uploaded once at init as static topology-like data) --
  // scaling them here is a plain host loop either way, same cost as
  // the cpuBackend path. Note this doesn't re-upload the (rare)
  // CUDA_Position_Restraint_Interaction's own static reference-position
  // copy -- a pre-existing gap, not introduced by this change (see
  // PERFORMANCE.md).
  if (scale_ref) {
    math::VArray & ref = conf.special().reference_positions;
    for (unsigned int i = 0; i < ref.size(); ++i)
      for (int j = 0; j < 3; ++j)
        ref(i)(j) *= mu(j);
  }

  // Own stream: lets configuration_view() insert cudaStreamWaitEvent()
  // against whichever GPU-native algorithm produced POS this step
  // (Leap_Frog_Position/constraints) instead of relying on legacy-
  // default-stream ordering.
  static cudaStream_t stream = 0;
  if (stream == 0) cudaStreamCreate(&stream);

  gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS, stream);
  const FPL9_TYPE mu9 = diag_to_fpl9(mu);
  gpu::launch_barostat_scale_positions(view.current().pos, num_atoms, mu9, stream);

  // Stay GPU-resident -- no eager sync-back. gpu_mirror_touches()
  // (berendsen_barostat.h) excludes POS from the framework's default
  // post-apply() invalidation, so this freshness survives; a later
  // CPU-only consumer or trajectory writeout publishes it lazily via
  // the existing flush_gpu_dirty()/needs_gpu_mirror_flush() machinery.
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_POS, stream);

  if (!sim.param().multicell.multicell && !math::boundary_check_cutoff(conf.current().box,
      sim.param().boundary.boundary, sim.param().pairlist.cutoff_long)) {
    io::messages.add("box is too small: not twice the cutoff!",
            "Berendsen_Barostat", io::message::error);
  }

  m_timer.stop();

  return 0;
}

template class algorithm::Berendsen_Barostat<util::gpuBackend>;
