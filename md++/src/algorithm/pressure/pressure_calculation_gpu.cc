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
 * @file pressure_calculation_gpu.cc
 * calculates the pressure -- GPU backend. Plain C++, no kernel syntax
 * (same split as berendsen_barostat_gpu.cc): the O(9) matrix combine is
 * trivially cheap regardless of backend and stays on the host, but
 * virial_tensor is read straight from the GPU mirror via
 * configuration_view() (the standard wait-on-fresh-producers protocol)
 * instead of the generic flush_gpu_dirty() path -- see
 * pressure_calculation.h's doc comment for why that distinction matters.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../configuration/state_properties.h"

#include "../../math/volume.h"
#include "../../math/transformation.h"

#include "../../gpu/cuda/manager/cuda_manager.h"

#include "pressure_calculation.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

#include "../../util/debug.h"

template<>
int algorithm::Pressure_Calculation<util::gpuBackend>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "Pressure calculation (GPU)");

  m_timer.start(sim);

  // Own stream: lets configuration_view() insert cudaStreamWaitEvent()
  // against whichever GPU-native algorithm produced VIRIAL this step
  // (bonded/nonbonded terms, CUDA_Molecular_Virial_Interaction,
  // CUDA_Lincs/CUDA_M_Shake/CUDA_Settle) instead of relying on legacy-
  // default-stream ordering.
  static cudaStream_t stream = 0;
  if (stream == 0) cudaStreamCreate(&stream);

  gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_VIRIAL, stream);
  cudaStreamSynchronize(stream);

  FPH9_TYPE vt_o;
  cudaMemcpy(&vt_o, view.old().virial_tensor, sizeof(FPH9_TYPE), cudaMemcpyDeviceToHost);

  math::Matrix virial;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      virial(i, j) = static_cast<double>(vt_o(i, j));

  // the virial is stored internally as just the outer product of positions and forces
  // so without the -0.5 prefactor.
  conf.old().pressure_tensor = topo.tot_cg_factor() *
          (conf.old().kinetic_energy_tensor + 0.5 * virial) *
          (2.0 / math::volume(conf.old().box, conf.boundary_type));

  m_timer.stop();

  return 0;

}

template class algorithm::Pressure_Calculation<util::gpuBackend>;
