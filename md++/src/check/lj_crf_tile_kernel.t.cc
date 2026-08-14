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
 * @file lj_crf_tile_kernel.t.cc
 * Standalone correctness test for gpu::lj_crf_tile_kernel (PLAN.md §10
 * step 8, TILE_PAIRLIST_DESIGN.md "step 8 onward"): kernel-only, not the
 * full CUDA_Nonbonded_Interaction wiring (that's step 9). Builds one real
 * 32x32 self-tile by hand from aladip's actual atoms/positions/charges/LJ
 * parameters (mirroring what classify_tiles_kernel would have produced:
 * mask bit set only for r < c pairs within [0.2nm, cutoff_long], see
 * compute_reference()'s doc comment) and compares the kernel's per-atom
 * forces and total LJ/CRF energies against
 * interaction::Nonbonded_Term::lj_crf_interaction (the same primitive
 * Standard_Pairlist_Algorithm's innerloop calls) evaluated directly on
 * the same pairs. Only USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../util/create_simulation.h"

#include "../math/periodicity.h"
#include "../interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "../interaction/nonbonded/interaction/nonbonded_term.h"

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/float3.h"
#include "gpu/cuda/memory/pairlist/tile.h"
#include "gpu/cuda/interaction/nonbonded/cuda_lj_params.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"
#include "gpu/cuda/interaction/nonbonded/kernels/lj_crf_tiles.h"

#include "check.h"

#ifdef XXMPI
  #include <mpi.h>
#endif

namespace {

  bool flush_messages(const char * label) {
    const io::message::severity_enum sev = io::messages.display(std::cerr);
    io::messages.clear();
    return sev >= io::message::error;
  }

  // Row/col block 0, col_from_b = false -- gpu::pack_block_index(0, 0,
  // false) (block_pairlist.h) encodes to exactly 0, but that header isn't
  // host-includable (it declares __global__ kernels taking
  // gpu::Periodicity<BOUNDARY>, which needs nvcc), so this is spelled out
  // instead of calling it.
  constexpr unsigned kSelfTileIndex = 0u;

  struct Reference {
    std::vector<math::Vec> force;
    double e_lj  = 0.0;
    double e_crf = 0.0;
  };

  /**
   * Builds one 32x32 self-tile over global atom indices [0, n), n =
   * min(32, num_atoms), *and* the matching CPU reference in the same
   * pass, via the exact primitive the CPU innerloop uses
   * (interaction::Nonbonded_Term::lj_crf_interaction) -- so both sides
   * are guaranteed to consider exactly the same pairs.
   *
   * mask bit (r, c) is set iff r < c (matching classify_tiles_kernel's
   * diagonal-tile convention: upper triangle only, no self-pair, no
   * double count) *and* the pair's distance falls in [min_dist,
   * cutoff_long] -- excluding pathologically close pairs (bonded
   * neighbors etc.) that a real tile would never contain anyway (they'd
   * be excluded upstream by classify_tiles_kernel) and whose enormous
   * r^-12 LJ repulsion amplifies float-vs-double rounding differences
   * far beyond what's meaningful for testing this kernel's arithmetic.
   */
  template <math::boundary_enum B>
  Reference compute_reference(topology::Topology & topo,
                               configuration::Configuration & conf,
                               simulation::Simulation & sim,
                               interaction::Nonbonded_Parameter & param,
                               unsigned n,
                               gpu::Interaction_Tile & tile) {
    interaction::Nonbonded_Term term;
    term.init(sim);
    math::Periodicity<B> periodicity(conf.current().box);

    const double min_dist2 = 0.2 * 0.2; // nm^2
    const double max_dist2 = sim.param().pairlist.cutoff_long *
                              sim.param().pairlist.cutoff_long;

    Reference ref;
    ref.force.resize(topo.num_atoms(), math::Vec(0.0, 0.0, 0.0));

    tile.index = kSelfTileIndex;
    for (unsigned r = 0; r < gpu::Interaction_Tile::MASK_SIZE; ++r) tile.mask[r] = 0u;

    for (unsigned i = 0; i < n; ++i) {
      for (unsigned j = i + 1; j < n; ++j) {
        math::Vec r;
        periodicity.nearest_image(conf.current().pos(i), conf.current().pos(j), r);
        const double dist2 = math::abs2(r);
        if (dist2 < min_dist2 || dist2 > max_dist2) continue;

        tile.mask[i] |= (1u << j);

        const interaction::lj_parameter_struct & lj = param.lj_parameter(topo.iac(i), topo.iac(j));
        const double q = topo.charge(i) * topo.charge(j);
        double f = 0.0, e_lj = 0.0, e_crf = 0.0;
        term.lj_crf_interaction(r, lj.c6, lj.c12, q, f, e_lj, e_crf);
        ref.force[i] += f * r;
        ref.force[j] -= f * r;
        ref.e_lj  += e_lj;
        ref.e_crf += e_crf;
      }
    }
    return ref;
  }

  int run_case(const std::string & stopo, const std::string & sconf,
               const std::string & sinput, math::boundary_enum force_boundary,
               const char * label, bool quiet) {
    util::simulation_struct s;
    io::In_Topology in_topo;
    in_topo.quiet = quiet;

    if (util::create_simulation(stopo, "", sconf, sinput, s, in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0) {
      std::cerr << label << ": creating simulation failed" << std::endl;
      return 1;
    }
    io::messages.display(std::cout);
    io::messages.clear();

    s.conf.boundary_type = force_boundary;
    s.sim.param().boundary.boundary = force_boundary;
    if (force_boundary == math::rectangular) {
      s.conf.current().box = math::Box(math::Vec(4.0, 0.0, 0.0),
                                        math::Vec(0.0, 4.0, 0.0),
                                        math::Vec(0.0, 0.0, 4.0));
    }

    const unsigned num_atoms = static_cast<unsigned>(s.topo.num_atoms());
    const unsigned n = std::min(gpu::Interaction_Tile::MASK_SIZE, num_atoms);

    // LJ parameter matrix, straight from the topology file (no forcefield
    // needed) -- interaction::Nonbonded_Parameter is the CPU-side
    // structure gpu::LJParams::init() reads from.
    std::vector<std::vector<interaction::lj_parameter_struct> > lj_matrix;
    in_topo.read_lj_parameter(lj_matrix, std::cout);
    interaction::Nonbonded_Parameter param;
    param.resize(static_cast<unsigned>(lj_matrix.size()));
    for (unsigned i = 0; i < lj_matrix.size(); ++i) {
      for (unsigned j = 0; j < lj_matrix[i].size(); ++j) {
        param.add_lj_parameter(i, j, lj_matrix[i][j]);
      }
    }
    if (flush_messages("lj parameter read")) {
      std::cerr << label << ": reading LJ parameters failed" << std::endl;
      return 1;
    }

    gpu::Interaction_Tile host_tile;
    Reference ref;
    if (force_boundary == math::vacuum) {
      ref = compute_reference<math::vacuum>(s.topo, s.conf, s.sim, param, n, host_tile);
    } else {
      ref = compute_reference<math::rectangular>(s.topo, s.conf, s.sim, param, n, host_tile);
    }

    // GPU side: same atoms/positions/charges/LJ params/tile, run through
    // gpu::lj_crf_tile_kernel via its runtime-boundary launch wrapper.
    gpu::LJParams gpu_lj;
    gpu_lj.init(param);

    interaction::Nonbonded_Term term;
    term.init(s.sim);
    gpu::NbSimParams nb;
    nb.four_pi_eps_i    = static_cast<FPL_TYPE>(math::four_pi_eps_i);
    nb.crf_2cut3i       = static_cast<FPL_TYPE>(term.crf_2cut3i());
    nb.crf_cut          = static_cast<FPL_TYPE>(term.crf_cut(0));
    nb.cutoff_short_sq  = static_cast<FPL_TYPE>(s.sim.param().pairlist.cutoff_short *
                                                 s.sim.param().pairlist.cutoff_short);
    nb.cutoff_long_sq   = static_cast<FPL_TYPE>(s.sim.param().pairlist.cutoff_long *
                                                 s.sim.param().pairlist.cutoff_long);

    math::CuVArray pos;
    pos.resize(num_atoms);
    gpu::cuvector<int> iac;
    iac.resize(num_atoms);
    gpu::cuvector<FPL_TYPE> charge;
    charge.resize(num_atoms);
    for (unsigned i = 0; i < num_atoms; ++i) {
      pos[i]    = static_cast<FPL3_TYPE>(s.conf.current().pos(i));
      iac[i]    = s.topo.iac(i);
      charge[i] = static_cast<FPL_TYPE>(s.topo.charge(i));
    }

    gpu::cuvector<unsigned> order;
    order.resize(n);
    for (unsigned i = 0; i < n; ++i) order[i] = i;

    gpu::TileVecT<gpu::Interaction_Tile> tiles;
    tiles.resize(1);
    tiles[0] = host_tile;

    gpu::cuvector<FPL3_TYPE> gpu_force;
    gpu_force.resize(num_atoms);
    for (unsigned i = 0; i < num_atoms; ++i) gpu_force[i] = make_FPL3(FPL_TYPE(0));

    gpu::cuvector<double> e_lj_total;
    e_lj_total.resize(1);
    e_lj_total[0] = 0.0;
    gpu::cuvector<double> e_crf_total;
    e_crf_total.resize(1);
    e_crf_total[0] = 0.0;

    gpu::launch_lj_crf_tiles(
        tiles.view(), order.data(), n, nullptr, 0u,
        pos.view(), iac.data(), charge.data(),
        gpu_lj.view(), nb, force_boundary, s.conf.current().box,
        gpu_force.data(), e_lj_total.data(), e_crf_total.data());
    cudaDeviceSynchronize();

    int errors = 0;
    // FPL_TYPE is float in the default (mixed-precision) build -- a
    // relative tolerance tighter than float's ~1e-7 precision would flag
    // ordinary float-vs-double rounding as a false failure.
    const double tol = 1e-4;

    if (std::abs(e_lj_total[0] - ref.e_lj) > tol * std::max(1.0, std::abs(ref.e_lj))) {
      std::cerr << label << ": e_lj mismatch: gpu=" << e_lj_total[0]
                << " cpu=" << ref.e_lj << std::endl;
      ++errors;
    }
    if (std::abs(e_crf_total[0] - ref.e_crf) > tol * std::max(1.0, std::abs(ref.e_crf))) {
      std::cerr << label << ": e_crf mismatch: gpu=" << e_crf_total[0]
                << " cpu=" << ref.e_crf << std::endl;
      ++errors;
    }
    for (unsigned i = 0; i < n; ++i) {
      const math::Vec gpu_f(gpu_force[i].x, gpu_force[i].y, gpu_force[i].z);
      const math::Vec diff = gpu_f - ref.force[i];
      const double scale = std::max(1.0, math::abs(ref.force[i]));
      if (math::abs(diff) > tol * scale) {
        std::cerr << label << ": force mismatch at atom " << i
                  << ": gpu=" << math::v2s(gpu_f)
                  << " cpu=" << math::v2s(ref.force[i]) << std::endl;
        ++errors;
      }
    }

    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK (" << n << " atoms, "
                << (n * (n - 1) / 2) << " pairs)" << std::endl;
    }
    return errors;
  }

} // namespace

int main(int argc, char* argv[]) {
#ifdef XXMPI
  MPI_Init(&argc, &argv);
#endif

  util::Known knowns;
  knowns << "verb";

  std::string usage = argv[0];
  usage += "\n\t[@verb   <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns, true)) {
    std::cerr << usage << std::endl;
    return 1;
  }
  util::parse_verbosity(args);
  const bool quiet = (args.count("verb") == -1);

  std::string stopo, sconf, sinput;
  GETFILEPATH(stopo, "aladip.topo", "src/check/data/");
  GETFILEPATH(sconf, "aladip.conf", "src/check/data/");
  GETFILEPATH(sinput, "aladip_unperturbed.in", "src/check/data/");

  int total = 0;
  total += run_case(stopo, sconf, sinput, math::vacuum, "vacuum", quiet);
  total += run_case(stopo, sconf, sinput, math::rectangular, "rectangular", quiet);

  return total;
}
