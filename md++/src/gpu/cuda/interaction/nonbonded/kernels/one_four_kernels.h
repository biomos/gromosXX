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
 * @file one_four_kernels.h
 * Host-callable entry point for 1,4-pair ("LJ exception") LJ+CRF
 * interactions -- the one nonbonded contribution missing from the GPU
 * tile pairlist entirely: `topo.all_exclusion(i)` (used for tile
 * classification masking, gpu::TopologyView::excl_ptr/excl_list) is
 * `topo.exclusion(i)` UNION `topo.one_four_pair(i)`
 * (`topology.cc`'s `update_all_exclusion()`), so 1,4 pairs are correctly
 * *excluded* from the regular tile-driven LJ/CRF sum on both CPU and
 * GPU -- but the CPU path adds them back in with their own scaled LJ
 * parameters (`cs6`/`cs12`) and Coulomb scaling factor via
 * `Nonbonded_Outerloop::one_four_outerloop` -> `one_four_interaction_
 * innerloop` (`nonbonded_innerloop.cc`), called unconditionally from
 * `Nonbonded_Set::calculate_interactions` (`nonbonded_set.cc`). Nothing
 * on the GPU path ever did this, so 1,4-pair-rich real topologies
 * (extended_test/ubiquitin: a real protein backbone, unlike aladip's
 * tiny test molecule) were silently missing this entire contribution --
 * found via a large, otherwise-unexplained solute-solute CRF mismatch
 * that survived ruling out the pairlist, precision, long-range caching,
 * and RF constants.
 *
 * Same host/device split as rf_excluded_kernels.h -- no `gpu/cuda/math/
 * periodicity.h` include, no `__global__` declaration, both private to
 * one_four_kernels.cu.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/interaction/nonbonded/cuda_lj_params.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"
#include "math/gmath.h"
#include "math/box.h"

namespace gpu {

  /**
   * @brief LJ (scaled cs6/cs12) + CRF (Coulomb-scaled) force/energy for
   * every 1,4 pair, added into `force`/`e_lj_total`/`e_crf_total`/
   * `virial_total` -- same accumulation convention (`+=`, per-[gi][gj]
   * energy-group bucket, flattened row-major 3x3 virial) as
   * `launch_lj_crf_tiles`/`launch_rf_excluded`, so all three are safe
   * to call back-to-back into the same buffers.
   *
   * `one_four_ptr`/`one_four_list`: CSR built from `topo.one_four_pair(i)`
   * (solute atoms only, j > i by construction -- same convention as
   * `rf_excl_ptr`/`rf_excl_list`, built the same way in
   * `CUDA_Pairlist_Algorithm_Impl::init()`).
   *
   * Exactly `interaction::Nonbonded_Term::lj_crf_interaction`'s formula
   * (`nonbonded_term.cc`) with `c6`/`c12` replaced by the pair's
   * `cs6`/`cs12` (`lj.get_scaled()`) and the real Coulomb term scaled by
   * `coulomb_scaling` (1.0 unless `param.amber.amber`, matching
   * `create_nonbonded.cc`'s `set_coulomb_scaling` call):
   *   e_crf = q*four_pi_eps_i*(coulomb_scaling/r - crf_2cut3i*r^2 - crf_cut)
   *
   * Does not synchronize -- same convention as the other nonbonded
   * kernels in this directory.
   */
  void launch_one_four(
      const int* one_four_ptr,
      const int* one_four_list,
      unsigned num_solute_atoms,
      math::CuVArray::View pos,
      const int* iac,
      const FPL_TYPE* charge,
      const unsigned* atom_energy_group,
      LJParamView lj,
      NbSimParams nb,
      FPL_TYPE coulomb_scaling,
      math::boundary_enum boundary,
      math::Box box,
      FPL3_TYPE* force,
      double* e_lj_total,
      double* e_crf_total,
      double* virial_total,
      cudaStream_t stream = nullptr);

}
