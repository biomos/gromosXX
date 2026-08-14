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
 * @file cuda_atomic_cutoff_toggle.t.cc
 * Regression test for the (formerly latent, now fixed) bug documented
 * in KNOWN_ISSUES.md as "CUDA context corruption after runtime
 * atomic_cutoff toggle": CUDA_Pairlist_Algorithm_Impl::prepare_cog()
 * sized m_cg_cog to num_solute_chargegroups() when it needed to cover
 * *all* chargegroups (prepare_cog_kernel writes cg_cog[cg_i] for every
 * chargegroup, and classify_tiles() reads it back for solvent
 * candidates too) -- an out-of-bounds write on every call whenever any
 * solvent chargegroups exist, found via `compute-sanitizer --tool
 * memcheck` while investigating why toggling atomic_cutoff seemed to
 * corrupt the CUDA context for the rest of the process (the OOB write
 * happens regardless of atomic_cutoff; it only sometimes clobbered
 * memory that mattered, which is why it looked toggle-specific).
 *
 * aladip_cuda.in has PERTURBATION on, so the real aladip_cuda ctest
 * never reaches CUDA_Pairlist_Algorithm::update()/prepare() and can't
 * exercise this at all. This drives the exact same sequence
 * (create_g96_forcefield -> ff->init() -> calculate_interactions() ->
 * toggle atomic_cutoff -> calculate_interactions() again, repeated for
 * several cycles) on a non-perturbed system (aladip_unperturbed.in)
 * with accelerator forced to gpu_cuda, then probes for corruption via
 * cudaGetLastError() and an unrelated fresh CUDA allocation. ctest
 * alone won't catch a silent out-of-bounds write that doesn't happen to
 * corrupt anything visible on a given run -- for that, run this binary
 * under `compute-sanitizer --tool memcheck` directly.
 */

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../interaction/interaction.h"
#include "../interaction/forcefield/forcefield.h"
#include "../interaction/forcefield/create_forcefield.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../util/create_simulation.h"

#include "check.h"

#include <cuda_runtime.h>

int main(int argc, char* argv[]) {
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

  util::simulation_struct s;
  io::In_Topology in_topo;
  in_topo.quiet = quiet;
  if (util::create_simulation(stopo, "", sconf, sinput, s, in_topo,
                               "", "", "", "", "", "", "", "", quiet) != 0) {
    std::cerr << "creating simulation failed" << std::endl;
    return 1;
  }
  io::messages.display(std::cout);
  io::messages.clear();

  s.sim.param().gpu.accelerator = simulation::gpu_cuda;
  s.sim.param().pcouple.virial = math::no_virial;

  interaction::Forcefield * ff = new interaction::Forcefield;
  if (interaction::create_g96_forcefield(*ff, s.topo, s.sim, in_topo, std::cout, quiet) != 0) {
    std::cerr << "creating forcefield failed" << std::endl;
    return 1;
  }
  io::messages.display(std::cout);
  io::messages.clear();

  if (ff->init(s.topo, s.conf, s.sim, std::cout, quiet) != 0) {
    std::cerr << "ff->init() failed" << std::endl;
    io::messages.display(std::cerr);
    return 1;
  }
  io::messages.display(std::cout);
  io::messages.clear();

  interaction::Interaction * nb = nullptr;
  for (auto it = ff->begin(); it != ff->end(); ++it) {
    std::cout << "interaction: " << (*it)->name << std::endl;
    if ((*it)->name == "NonBonded") nb = *it;
  }
  if (!nb) {
    std::cerr << "no NonBonded interaction found" << std::endl;
    return 1;
  }

  auto report_cuda_error = [](const char * label) {
    cudaError_t err = cudaGetLastError();
    std::cout << label << ": cudaGetLastError() = " << cudaGetErrorString(err)
              << " (" << static_cast<int>(err) << ")" << std::endl;
    return err;
  };

  std::cout << "--- step 1: calculate_interactions(), atomic_cutoff=false ---" << std::endl;
  s.conf.current().force = 0.0;
  s.conf.current().energies.zero();
  int ret1 = nb->calculate_interactions(s.topo, s.conf, s.sim);
  std::cout << "ret1=" << ret1 << std::endl;
  io::messages.display(std::cout);
  io::messages.clear();
  report_cuda_error("after step 1");

  cudaError_t err2 = cudaSuccess, err3 = cudaSuccess;
  for (int cycle = 0; cycle < 10; ++cycle) {
    std::cout << "--- cycle " << cycle << ": toggle atomic_cutoff=true, calculate_interactions() ---" << std::endl;
    s.sim.param().pairlist.atomic_cutoff = true;
    s.conf.current().force = 0.0;
    s.conf.current().energies.zero();
    int ret2 = nb->calculate_interactions(s.topo, s.conf, s.sim);
    std::cout << "ret2=" << ret2 << std::endl;
    io::messages.display(std::cout);
    io::messages.clear();
    err2 = report_cuda_error("after toggle true");
    if (err2 != cudaSuccess) break;

    std::cout << "--- cycle " << cycle << ": toggle atomic_cutoff=false, calculate_interactions() ---" << std::endl;
    s.sim.param().pairlist.atomic_cutoff = false;
    s.conf.current().force = 0.0;
    s.conf.current().energies.zero();
    int ret3 = nb->calculate_interactions(s.topo, s.conf, s.sim);
    std::cout << "ret3=" << ret3 << std::endl;
    io::messages.display(std::cout);
    io::messages.clear();
    err3 = report_cuda_error("after toggle false");
    if (err3 != cudaSuccess) break;
  }

  std::cout << "--- step 4: unrelated fresh CUDA op (cudaMalloc/cudaFree) ---" << std::endl;
  void * p = nullptr;
  cudaError_t merr = cudaMalloc(&p, 1024);
  std::cout << "cudaMalloc: " << cudaGetErrorString(merr) << std::endl;
  if (merr == cudaSuccess) cudaFree(p);
  report_cuda_error("after step 4 (unrelated cudaMalloc)");

  std::cout << "--- step 5: cudaDeviceSynchronize probe ---" << std::endl;
  cudaError_t serr = cudaDeviceSynchronize();
  std::cout << "cudaDeviceSynchronize: " << cudaGetErrorString(serr) << std::endl;

  delete ff;

  const bool corrupted = (err2 != cudaSuccess) || (err3 != cudaSuccess) ||
                          (merr != cudaSuccess) || (serr != cudaSuccess);
  std::cout << (corrupted ? "REPRO: corruption observed" : "no corruption observed") << std::endl;
  return corrupted ? 1 : 0;
}
