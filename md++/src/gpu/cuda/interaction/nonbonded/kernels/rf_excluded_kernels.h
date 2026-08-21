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
 * @file rf_excluded_kernels.h
 * Host-callable entry point for the reaction-field-for-excluded-pairs
 * kernel (`param.nonbonded.rf_excluded`), the one nonbonded contribution
 * `gpu::launch_lj_crf_tiles`'s tile pairlist deliberately never sees --
 * excluded pairs are masked out of every tile during classification
 * (`classify_tiles_kernel`, block_pairlist.cu), so this term needs its
 * own kernel that walks the exclusion list directly, exactly mirroring
 * `Nonbonded_Outerloop::RF_excluded_outerloop`'s two CPU innerloops
 * (`nonbonded_innerloop.cc`):
 *
 *   - `RF_excluded_interaction_innerloop` (solute): for every solute atom
 *     i, a self-term (r=0, energy only) plus one real pair term per
 *     entry in `topo.exclusion(i)`.
 *   - `RF_solvent_interaction_innerloop` (solvent): for every solvent
 *     chargegroup (== one rigid molecule), the distance-dependent energy
 *     term for every pair of atoms *within* that chargegroup -- no self
 *     term (its distance-independent part is left out deliberately, see
 *     the CPU comment: "should add up to zero"), no force (rigid
 *     solvent).
 *
 * Same host/device split as lj_crf_tiles.h: no `gpu/cuda/math/
 * periodicity.h` include here, no `__global__` declaration -- both
 * kernels are private to rf_excluded_kernels.cu, which instantiates them
 * directly from launch_rf_excluded()'s boundary switch.
 *
 * Unlike the tile kernel, this is not called through the pairlist at
 * all -- it's independent of pairlist_update/skip_step, called every
 * step directly from CUDA_Nonbonded_Interaction::calculate_interactions,
 * matching nonbonded_set.cc's unconditional (non-twin-range) call.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"
#include "math/gmath.h"
#include "math/box.h"

namespace gpu {

  /**
   * @brief Reaction-field self-term + excluded-pair energy/force
   * (solute) and rigid-solvent excluded-pair energy (solvent), added
   * into `force`/`e_crf_total`/`virial_total` -- same accumulation
   * convention (`+=`, per-[gi][gj] energy-group bucket, flattened
   * row-major 3x3 virial) as `launch_lj_crf_tiles`, so the two are safe
   * to call back-to-back into the same buffers without an intervening
   * zero.
   *
   * `rf_excl_ptr`/`rf_excl_list`: CSR built from `topo.exclusion(i)`
   * (plain exclusions only, solute atoms only, size
   * num_solute_atoms+1/num_exclusion_entries) -- deliberately NOT
   * `gpu::TopologyView::excl_ptr`/`excl_list` (which is built from
   * `topo.all_exclusion(i)` = `exclusion(i)` UNION `one_four_pair(i)`,
   * `topology.cc`'s `update_all_exclusion()`): 1-4 pairs get their own
   * separate `lj_exception` treatment, never the RF correction, and
   * reusing the tile-classification CSR here would silently apply it to
   * them too. `CUDA_Pairlist_Algorithm_Impl::init()` builds this
   * dedicated CSR once (`m_rf_excl_ptr`/`m_rf_excl_list`), the same way
   * it builds `m_iac`/`m_charge`/`m_atom_energy_group`.
   *
   * `chargegroup`/`num_solute_chargegroups`/`num_chargegroups`: the CSR
   * chargegroup-offset array `gpu::TopologyView` already carries
   * (solvent chargegroups occupy
   * [num_solute_chargegroups, num_chargegroups)) -- fine to reuse as-is,
   * unrelated to the exclusion-vs-1-4 distinction above.
   * `charge`/`atom_energy_group` are the same `FPL_TYPE`/`unsigned`
   * arrays `CUDA_Pairlist_Algorithm_Impl` already builds and passes to
   * `launch_lj_crf_tiles` (`m_charge`/`m_atom_energy_group`).
   *
   * Does not synchronize -- same convention as launch_lj_crf_tiles.
   */
  void launch_rf_excluded(
      const int* rf_excl_ptr,
      const int* rf_excl_list,
      unsigned num_solute_atoms,
      const int* chargegroup,
      unsigned num_solute_chargegroups,
      unsigned num_chargegroups,
      math::CuVArray::View pos,
      const FPL_TYPE* charge,
      const unsigned* atom_energy_group,
      NbSimParams nb,
      math::boundary_enum boundary,
      math::Box box,
      FPH3_TYPE* force,
      double* e_crf_total,
      double* virial_total,
      cudaStream_t stream = nullptr);

}
