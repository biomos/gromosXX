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

#include "../../stdheader.h"

#include <cstdlib>
#include <fstream>
#include <memory>
#include <sys/stat.h>

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../algorithm/algorithm/algorithm_sequence.h"
#include "../../io/argument.h"
#include "../../io/message.h"
#include "../../io/read_input.h"
#include "../../io/configuration/out_configuration.h"

#include "ubiquitin_runner.h"

namespace {
  bool flush_messages(std::ostream & err, const char * label) {
    const io::message::severity_enum sev = io::messages.display(err);
    io::messages.clear();
    if (sev >= io::message::error) {
      err << "ubiquitin runner: " << label << " logged an error-or-worse message"
          << std::endl;
      return true;
    }
    return false;
  }
}

int check_ubiquitin::run_ubiquitin(const std::string & imd_path, unsigned num_steps,
                                    std::vector<check_ubiquitin::StepResult> & out,
                                    std::ostream & err,
                                    const check_ubiquitin::OutputPaths * out_paths) {
  std::string stopo, sconf;
  {
    struct stat buffer;
    stopo = TOP_SOURCE_DIR "/src/check/extended_test/ubiquitin_54a8.top";
    if (stat(stopo.c_str(), &buffer) != 0)
      stopo = "../../" TOP_SOURCE_DIR "/src/check/extended_test/ubiquitin_54a8.top";
    sconf = TOP_SOURCE_DIR "/src/check/extended_test/ubiquitin.cnf";
    if (stat(sconf.c_str(), &buffer) != 0)
      sconf = "../../" TOP_SOURCE_DIR "/src/check/extended_test/ubiquitin.cnf";
  }

  io::Argument args;
  args.insert(std::make_pair(std::string("topo"), stopo));
  args.insert(std::make_pair(std::string("conf"), sconf));
  args.insert(std::make_pair(std::string("input"), imd_path));
  if (out_paths) {
    args.insert(std::make_pair(std::string("fin"), out_paths->fin));
    args.insert(std::make_pair(std::string("trc"), out_paths->trc));
    args.insert(std::make_pair(std::string("tre"), out_paths->tre));
  }

  topology::Topology topo;
  configuration::Configuration conf;
  simulation::Simulation sim;
  algorithm::Algorithm_Sequence md_seq;

  std::ostringstream discard;
  if (io::read_input(args, topo, conf, sim, md_seq, discard, true)) {
    flush_messages(err, "read_input");
    err << "ubiquitin runner: read_input failed" << std::endl;
    return 1;
  }
  if (flush_messages(err, "read_input")) return 1;

  // Real io::Out_Configuration, driven exactly like program/md.cc's own
  // main loop (traj.write() before each run(), traj.print() after,
  // traj.write(..., io::final)/traj.print_final() once at the end) --
  // see ubiquitin_runner.h's doc comment for why this matters as its
  // own, independently-checked code path.
  std::unique_ptr<io::Out_Configuration> traj;
  if (out_paths) {
    traj.reset(new io::Out_Configuration(GROMOSXX "\n"));
    traj->title(GROMOSXX "\n" + sim.param().title);
    traj->init(args, sim.param());
    if (flush_messages(err, "Out_Configuration::init")) return 1;
  }

  if (md_seq.init(topo, conf, sim, discard, true)) {
    flush_messages(err, "Algorithm_Sequence::init");
    err << "ubiquitin runner: Algorithm_Sequence::init failed" << std::endl;
    return 1;
  }
  if (flush_messages(err, "Algorithm_Sequence::init")) return 1;

  out.clear();
  out.reserve(num_steps);

  for (unsigned step = 0; step < num_steps; ++step) {
    if (traj) traj->write(conf, topo, sim, io::reduced);

    if (md_seq.run(topo, conf, sim)) {
      flush_messages(err, "Algorithm_Sequence::run");
      err << "ubiquitin runner: run failed at step " << step << std::endl;
      return 1;
    }
    if (flush_messages(err, "Algorithm_Sequence::run")) return 1;

    if (traj) traj->print(topo, conf, sim);

    // Forcefield::calculate_interactions writes into conf.current()
    // (see e.g. cuda_pairlist_algorithm_impl.cu's conf.current().
    // energies.lj_energy[...] += ... accumulation), but some algorithm
    // in the sequence (leap-frog's own old/current rotation) has
    // already moved that into conf.old() by the time this step's
    // run() call returns -- confirmed directly against
    // io::Out_Configuration::print(), which reads conf.old().energies
    // (out_configuration.cc), not conf.current(). Reading
    // conf.current() here instead left every in-memory StepResult
    // silently all-zero, invisible until this same value was cross-
    // checked against a real .tre file's contents (see this file's
    // own doc comment for why that cross-check exists at all).
    StepResult r;
    r.step = step;
    const configuration::Energy & e = conf.old().energies;
    r.total             = e.total;
    r.kinetic_total      = e.kinetic_total;
    r.potential_total    = e.potential_total;
    r.bonded_total        = e.bonded_total;
    r.nonbonded_total    = e.nonbonded_total;
    r.lj_total            = e.lj_total;
    r.crf_total           = e.crf_total;
    r.constraints_total  = e.constraints_total;

    r.baths.clear();
    for (unsigned b = 0; b < sim.multibath().size(); ++b) {
      BathEnergy be;
      be.ekin = sim.multibath()[b].ekin;
      r.baths.push_back(be);
    }
    out.push_back(r);

    sim.steps() = sim.steps() + sim.param().analyze.stride;
    sim.time() = sim.param().step.t0 + sim.steps() * sim.time_step_size();
  }

  if (traj) {
    traj->write(conf, topo, sim, io::final);
    traj->print_final(topo, conf, sim);
  }

  return 0;
}

void check_ubiquitin::write_reference(const std::string & path,
                                       const std::vector<check_ubiquitin::StepResult> & steps) {
  std::ofstream f(path);
  f.precision(17);
  f << std::scientific;
  for (const StepResult & r : steps) {
    f << r.step << ' ' << r.total << ' ' << r.kinetic_total << ' '
      << r.potential_total << ' ' << r.bonded_total << ' ' << r.nonbonded_total
      << ' ' << r.lj_total << ' ' << r.crf_total << ' ' << r.constraints_total
      << ' ' << r.baths.size();
    for (const BathEnergy & b : r.baths) f << ' ' << b.ekin;
    f << '\n';
  }
}

bool check_ubiquitin::read_reference(const std::string & path,
                                      std::vector<check_ubiquitin::StepResult> & out,
                                      std::ostream & err) {
  std::ifstream f(path);
  if (!f) {
    err << "ubiquitin runner: could not open reference file " << path << std::endl;
    return false;
  }
  out.clear();
  StepResult r;
  unsigned num_baths;
  while (f >> r.step >> r.total >> r.kinetic_total >> r.potential_total
           >> r.bonded_total >> r.nonbonded_total >> r.lj_total >> r.crf_total
           >> r.constraints_total >> num_baths) {
    r.baths.assign(num_baths, BathEnergy{0.0});
    for (unsigned b = 0; b < num_baths; ++b) f >> r.baths[b].ekin;
    out.push_back(r);
  }
  if (out.empty()) {
    err << "ubiquitin runner: reference file " << path << " is empty or malformed"
        << std::endl;
    return false;
  }
  return true;
}

std::string check_ubiquitin::make_temp_dir(std::ostream & err) {
  std::string tmpl = "/tmp/ubiquitin_test_XXXXXX";
  std::vector<char> buf(tmpl.begin(), tmpl.end());
  buf.push_back('\0');
  if (mkdtemp(buf.data()) == nullptr) {
    err << "ubiquitin runner: mkdtemp failed" << std::endl;
    return "";
  }
  return std::string(buf.data());
}
