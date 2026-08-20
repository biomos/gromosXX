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
 * @file temperature_gpu.t.cc
 * Validates Temperature_Calculation<gpuBackend> and
 * Berendsen_Thermostat<gpuBackend> (PLAN.md §10 step 13) against the
 * CPU reference, on a synthetic multibath setup that deliberately
 * exercises BOTH of Thermostat::scale()'s CPU-side cases in one run:
 * a "separately coupled" temperature group (com_bath != ir_bath) and a
 * "jointly coupled" one (com_bath == ir_bath) -- aladip's own MULTIBATH
 * block only ever uses com_bath == ir_bath, so the temperature groups
 * and bath assignment are overridden here the same way
 * cuda_nonbonded_interaction.t.cc overrides energy groups. Only
 * USE_CUDA builds run this.
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
#include "../algorithm/temperature/berendsen_thermostat.h"

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

  /**
   * Overrides topo.temperature_groups()/sim.multibath() with a
   * synthetic 2-temperature-group, 3-bath setup: group 0 (first half
   * of the atoms) is separately coupled (com_bath=0, ir_bath=1); group
   * 1 (second half) is jointly coupled (com_bath=2 == ir_bath=2).
   * Applied identically to both the CPU and GPU simulation_structs.
   */
  void setup_synthetic_multibath(util::simulation_struct & s) {
    const unsigned num_atoms = static_cast<unsigned>(s.topo.num_atoms());
    const unsigned mid = num_atoms / 2;

    // topo.temperature_groups() is NOT the same convention as
    // energy_groups(): it's num_groups+1 entries with a leading 0 and
    // EXCLUSIVE end boundaries (confirmed directly from in_topology.cc's
    // TEMPERATUREGROUPS parsing: "if size()==0 { push_back(0);
    // push_back(num_solute_atoms); }", num_solute_temperature_groups() =
    // size()-1). Group g spans atoms [tg[g], tg[g+1]).
    s.topo.temperature_groups().assign({0u, mid, num_atoms});

    s.sim.multibath().clear();
    s.sim.multibath().add_bath(298.15, 0.01, 100.0, 30.0, 70.0);   // bath 0
    s.sim.multibath().add_bath(305.0,  0.01, 100.0, 30.0, 70.0);   // bath 1
    s.sim.multibath().add_bath(298.15, 0.01, 200.0, 60.0, 140.0);  // bath 2

    // bath_index_struct's last_atom IS 0-indexed inclusive (confirmed
    // from in_parameter.cc's MULTIBATH parsing: "add_bath_index(last-1,
    // ..., com_bath-1, ir_bath-1)"), matching energy_groups()'s
    // convention -- unlike temperature_groups() above. last_temperature_
    // group is also 0-indexed (multibath.cc: "last_temperature_group =
    // tg - 1"); each range below covers exactly one temperature group.
    s.sim.multibath().bath_index().clear();
    s.sim.multibath().bath_index().push_back(
        simulation::bath_index_struct(mid - 1, /*last_temp_group=*/0u, /*com=*/0u, /*ir=*/1u));
    s.sim.multibath().bath_index().push_back(
        simulation::bath_index_struct(num_atoms - 1, /*last_temp_group=*/1u, /*com=*/2u, /*ir=*/2u));

    // The file's own MULTIBATH block only ever defines 2 baths --
    // conf.*().energies.kinetic_energy/com_kinetic_energy/
    // ir_kinetic_energy were sized for that at simulation-construction
    // time. Resize for the 3-bath synthetic setup above (keeping the
    // existing energy-group count), or writes to bath index 2 are
    // out-of-bounds.
    const unsigned num_energy_groups =
        static_cast<unsigned>(s.conf.current().energies.bond_energy.size());
    s.conf.current().energies.resize(num_energy_groups, 3);
    s.conf.old().energies.resize(num_energy_groups, 3);
  }

  int compare_energies(const util::simulation_struct & cpu_s,
                        const util::simulation_struct & gpu_s,
                        const char * label, double tol) {
    int errors = 0;
    const unsigned num_baths = static_cast<unsigned>(cpu_s.sim.multibath().size());
    for (unsigned b = 0; b < num_baths; ++b) {
      const double cpu_ekin = cpu_s.sim.multibath()[b].ekin;
      const double gpu_ekin = gpu_s.sim.multibath()[b].ekin;
      if (std::abs(cpu_ekin - gpu_ekin) > tol * std::max(1.0, std::abs(cpu_ekin))) {
        std::cerr << label << ": bath " << b << " ekin mismatch: cpu=" << cpu_ekin
                  << " gpu=" << gpu_ekin << std::endl;
        ++errors;
      }
      const double cpu_com = cpu_s.conf.old().energies.com_kinetic_energy[b];
      const double gpu_com = gpu_s.conf.old().energies.com_kinetic_energy[b];
      if (std::abs(cpu_com - gpu_com) > tol * std::max(1.0, std::abs(cpu_com))) {
        std::cerr << label << ": bath " << b << " com_kinetic_energy mismatch: cpu="
                  << cpu_com << " gpu=" << gpu_com << std::endl;
        ++errors;
      }
      const double cpu_ir = cpu_s.conf.old().energies.ir_kinetic_energy[b];
      const double gpu_ir = gpu_s.conf.old().energies.ir_kinetic_energy[b];
      if (std::abs(cpu_ir - gpu_ir) > tol * std::max(1.0, std::abs(cpu_ir))) {
        std::cerr << label << ": bath " << b << " ir_kinetic_energy mismatch: cpu="
                  << cpu_ir << " gpu=" << gpu_ir << std::endl;
        ++errors;
      }
    }
    return errors;
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

  int run_case(const std::string & stopo, const std::string & sconf,
               const std::string & sinput, const char * label, bool quiet) {
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

    const double tol = 1e-4;
    int errors = 0;

    // --- Temperature_Calculation ---
    algorithm::Temperature_Calculation<util::cpuBackend> cpu_tcalc;
    algorithm::Temperature_Calculation<util::gpuBackend> gpu_tcalc;

    if (cpu_tcalc.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_tcalc.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("temperature_calculation apply");
      std::cerr << label << ": Temperature_Calculation::apply() failed" << std::endl;
      return 1;
    }
    // gpu_tcalc's apply() only launches the reduction kernels -- the
    // per-bath bookkeeping (bath.ekin, conf.old().energies.*) is
    // deferred to finalize_gpu_step(), normally resolved by
    // Algorithm_Sequence::run() (see Algorithm::finalize_gpu_step()'s
    // doc comment). This standalone test drives apply() directly, so
    // it must call it explicitly before reading the result below.
    gpu_tcalc.finalize_gpu_step(gpu_s.topo, gpu_s.conf, gpu_s.sim);
    flush_messages("temperature_calculation apply");

    errors += compare_energies(cpu_s, gpu_s, (std::string(label) + " (temperature)").c_str(), tol);

    // --- Berendsen_Thermostat (uses the bath.ekin Temperature_Calculation
    // just set) ---
    algorithm::Berendsen_Thermostat<util::cpuBackend> cpu_thermo;
    algorithm::Berendsen_Thermostat<util::gpuBackend> gpu_thermo;

    if (cpu_thermo.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_thermo.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("berendsen_thermostat apply");
      std::cerr << label << ": Berendsen_Thermostat::apply() failed" << std::endl;
      return 1;
    }
    flush_messages("berendsen_thermostat apply");

    // Berendsen_Thermostat<gpuBackend> deliberately leaves the scaled
    // velocity resident on the GPU mirror (mark_gpu_dirty(), no
    // sync-back) -- in the real sequence, Leap_Frog_Position<gpuBackend>
    // running right after it is what eventually flushes to the CPU.
    // This test drives the thermostat standalone, so it must do that
    // sync itself before comparing conf.current().vel() directly.
    gpu_s.sim.cuda().sync_configuration_from_device(gpu_s.conf);

    errors += compare_velocities(cpu_s, gpu_s, (std::string(label) + " (thermostat)").c_str(), tol);

    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK (3 baths, com!=ir + com==ir coupling both exercised)"
                << std::endl;
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

  return run_case(stopo, sconf, sinput, "temperature_gpu", quiet);
}
