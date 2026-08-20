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
 * @file nosehoover_gpu.t.cc
 * Validates NoseHoover_Thermostat<gpuBackend> (PLAN.md §10 step 20)
 * against the CPU reference. Same synthetic 2-temperature-group/
 * 3-bath setup as temperature_gpu.t.cc (deliberately exercising both
 * of Thermostat::scale()'s cases, since NoseHoover_Thermostat<
 * gpuBackend> reuses that exact kernel) -- run twice, once with plain
 * Nose-Hoover (multibath.algorithm == 1, calc_scaling()) and once with
 * a Nose-Hoover chain (multibath.algorithm == 3, calc_chain_scaling()),
 * since those are the two backend-agnostic scalar-math paths apply()
 * dispatches between; both must feed the same GPU velocity-scaling
 * kernel correctly. Only USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../simulation/multibath.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/temperature/temperature_calculation.h"
#include "../algorithm/temperature/nosehoover_thermostat.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../util/create_simulation.h"

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

  // Identical to temperature_gpu.t.cc's setup_synthetic_multibath():
  // group 0 (first half of the atoms) is separately coupled
  // (com_bath=0, ir_bath=1); group 1 (second half) is jointly coupled
  // (com_bath=2 == ir_bath=2).
  void setup_synthetic_multibath(util::simulation_struct & s) {
    const unsigned num_atoms = static_cast<unsigned>(s.topo.num_atoms());
    const unsigned mid = num_atoms / 2;

    s.topo.temperature_groups().assign({0u, mid, num_atoms});

    s.sim.multibath().clear();
    s.sim.multibath().add_bath(298.15, 0.01, 100.0, 30.0, 70.0);   // bath 0
    s.sim.multibath().add_bath(305.0,  0.01, 100.0, 30.0, 70.0);   // bath 1
    s.sim.multibath().add_bath(298.15, 0.01, 200.0, 60.0, 140.0);  // bath 2

    s.sim.multibath().bath_index().clear();
    s.sim.multibath().bath_index().push_back(
        simulation::bath_index_struct(mid - 1, /*last_temp_group=*/0u, /*com=*/0u, /*ir=*/1u));
    s.sim.multibath().bath_index().push_back(
        simulation::bath_index_struct(num_atoms - 1, /*last_temp_group=*/1u, /*com=*/2u, /*ir=*/2u));

    const unsigned num_energy_groups =
        static_cast<unsigned>(s.conf.current().energies.bond_energy.size());
    s.conf.current().energies.resize(num_energy_groups, 3);
    s.conf.old().energies.resize(num_energy_groups, 3);
  }

  int compare_velocities(const util::simulation_struct & cpu_s,
                          const util::simulation_struct & gpu_s,
                          const char * label, double tol) {
    int errors = 0;
    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());
    for (unsigned i = 0; i < num_atoms; ++i) {
      const math::Vec diff = cpu_s.conf.current().vel(i) - gpu_s.conf.current().vel(i);
      const double scale = std::max(1.0, math::abs(cpu_s.conf.current().vel(i)));
      if (math::abs(diff) > tol * scale) {
        std::cerr << label << ": vel mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().vel(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().vel(i)) << std::endl;
        ++errors;
      }
    }
    return errors;
  }

  int run_one(const std::string & stopo, const std::string & sconf,
              const std::string & sinput, int multibath_algorithm,
              const char * label, bool quiet) {
    util::simulation_struct cpu_s, gpu_s;
    io::In_Topology cpu_in_topo, gpu_in_topo;
    cpu_in_topo.quiet = quiet;
    gpu_in_topo.quiet = quiet;

    if (util::create_simulation(stopo, "", sconf, sinput, cpu_s, cpu_in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0 ||
        util::create_simulation(stopo, "", sconf, sinput, gpu_s, gpu_in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0) {
      std::cerr << label << ": creating simulation failed" << std::endl;
      return 1;
    }
    io::messages.display(std::cout);
    io::messages.clear();

    setup_synthetic_multibath(cpu_s);
    setup_synthetic_multibath(gpu_s);

    cpu_s.sim.param().multibath.algorithm = multibath_algorithm;
    gpu_s.sim.param().multibath.algorithm = multibath_algorithm;

    const double tol = 1e-4;
    int errors = 0;

    // --- Temperature_Calculation (sets bath.ekin, needed by
    // calc_scaling()/calc_chain_scaling()) ---
    algorithm::Temperature_Calculation<util::cpuBackend> cpu_tcalc;
    algorithm::Temperature_Calculation<util::gpuBackend> gpu_tcalc;

    if (cpu_tcalc.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_tcalc.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("temperature_calculation apply");
      std::cerr << label << ": Temperature_Calculation::apply() failed" << std::endl;
      return 1;
    }
    // gpu_tcalc's apply() only launches the reduction kernels -- see
    // temperature_gpu.t.cc's identical comment for why this standalone
    // test needs its own explicit finalize_gpu_step() call.
    gpu_tcalc.finalize_gpu_step(gpu_s.topo, gpu_s.conf, gpu_s.sim);
    flush_messages("temperature_calculation apply");

    // --- NoseHoover_Thermostat ---
    algorithm::NoseHoover_Thermostat<util::cpuBackend> cpu_nh;
    algorithm::NoseHoover_Thermostat<util::gpuBackend> gpu_nh;

    if (cpu_nh.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_nh.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("nosehoover_thermostat init");
      std::cerr << label << ": NoseHoover_Thermostat::init() failed" << std::endl;
      return 1;
    }
    flush_messages("nosehoover_thermostat init");

    if (cpu_nh.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_nh.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("nosehoover_thermostat apply");
      std::cerr << label << ": NoseHoover_Thermostat::apply() failed" << std::endl;
      return 1;
    }
    flush_messages("nosehoover_thermostat apply");

    // NoseHoover_Thermostat<gpuBackend> deliberately leaves the scaled
    // velocity resident on the GPU mirror (see berendsen_thermostat_gpu.
    // cc's identical reasoning) -- flush it explicitly for this
    // standalone comparison (see temperature_gpu.t.cc's identical note).
    gpu_s.sim.cuda().sync_configuration_from_device(gpu_s.conf);

    errors += compare_velocities(cpu_s, gpu_s, label, tol);

    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK" << std::endl;
    }
    return errors ? 1 : 0;
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

  int errors = 0;
  errors += run_one(stopo, sconf, sinput, 1, "nosehoover_gpu (plain)", quiet);
  errors += run_one(stopo, sconf, sinput, 3, "nosehoover_gpu (chain, 3 instances)", quiet);
  return errors ? 1 : 0;
}
