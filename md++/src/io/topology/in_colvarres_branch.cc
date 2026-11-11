/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres_branch.cc
 * implements methods of In_Topology for COLVAR geometry specification files.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"

#include "in_colvarres.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

#include "in_colvarres_helpers.h"

using namespace io::colvar_input;

namespace io {
namespace colvar_input {
void read_colvar_branches(const std::map<std::string, std::vector<std::string> > &blocks,
    simulation::Simulation &sim, std::set<std::string> &read, std::ostream &os, bool quiet)
{
  for (bool perturbed : {false, true}) {
    const std::string name = perturbed ? "PERTDIHEDRALBRANCH" : "DIHEDRALBRANCH";
    const auto block = blocks.find(name);
    if (block == blocks.end() || block->second.size() <= 2) continue;
    read.insert(name);
    auto &specs = perturbed ? sim.param().colvarres.pert_bias_specs : sim.param().colvarres.bias_specs;
    std::set<unsigned int> used;
    if (!quiet) os << name << "\n# INDX DELTA (degrees)\n";
    for (auto line = block->second.begin()+1; line != block->second.end()-1; ++line) {
      if (!line_is_data(*line)) continue;
      long long index = 0;
      double delta = 0;
      std::string extra;
      std::istringstream input(line->substr(0, line->find('#')));
      if (!(input >> index >> delta) || (input >> extra) || index <= 0
          || index > std::numeric_limits<unsigned int>::max()
          || !std::isfinite(delta)) {
        io::messages.add(name + ": expected positive INDX and finite DELTA in degrees.",
                         "In_Colvarres", io::message::error);
        continue;
      }
      if (!used.insert(static_cast<unsigned int>(index)).second) {
        io::messages.add(name + ": duplicate INDX.", "In_Colvarres", io::message::error);
        continue;
      }
      // Reduce in degrees before conversion to avoid overflow and loss of
      // turn information for very large but finite input offsets.
      double canonical = std::fmod(delta, 360.0);
      if (canonical < 0.0) canonical += 360.0;
      if (canonical == 0.0 && delta > 0.0) canonical = 360.0;
      bool found = false;
      for (auto &spec : specs) {
        if (spec.index != static_cast<unsigned int>(index)) continue;
        found = true;
        if (spec.type != "DIHEDRAL" || spec.rah != 0) {
          io::messages.add(name + ": explicit branches require a full DIHEDRAL restraint (RAH 0).",
                           "In_Colvarres", io::message::error);
          continue;
        }
        spec.explicit_periodic_branch = true;
        spec.branch_delta = canonical * math::Pi / 180;
      }
      if (!found) io::messages.add(name + ": INDX has no matching restraint.", "In_Colvarres", io::message::error);
      if (!quiet && found) os << index << ' ' << delta << '\n';
    }
    if (!quiet) os << "END\n";
  }
}
} // namespace colvar_input
} // namespace io
