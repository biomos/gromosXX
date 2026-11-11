/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres_eds.cc
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

void read_eds_colvar_blocks(const std::map<std::string, std::vector<std::string> > &blocks,
    topology::Topology &topo, simulation::Simulation &sim, std::set<std::string> &block_read)
{
  const auto data = [&](const std::string &name) {
    std::vector<std::string> result;
    const auto found = blocks.find(name);
    if (found == blocks.end() || found->second.size() <= 2) return result;
    block_read.insert(name);
    for (auto it = found->second.begin()+1; it != found->second.end()-1; ++it)
      if (line_is_data(*it)) result.push_back(it->substr(0,it->find('#')));
    return result;
  };
  const auto error = [](const std::string &message) {
    io::messages.add(message, "In_Colvarres", io::message::error);
  };
  const auto rows = data("EDSRESTRAINTS");
  const auto geometry = data("EDSDISTANCERESSPEC");
  if (rows.empty() && geometry.empty()) return;
  if (!sim.param().eds.eds || sim.param().eds.numstates < 2) {
    error("EDSRESTRAINTS requires enabled EDS/AEDS with at least two states.");
    return;
  }
  if (sim.param().eds.form < simulation::single_s || sim.param().eds.form > simulation::aeds_advanced_search2) {
    error("EDSRESTRAINTS requires a supported EDS or AEDS Hamiltonian.");
    return;
  }
  if (rows.empty() || geometry.empty()) {
    error("EDSRESTRAINTS and EDSDISTANCERESSPEC must both be supplied.");
    return;
  }
  unsigned int count = 0;
  std::string extra;
  std::istringstream first(rows.front());
  if (!(first >> count) || (first >> extra) || count == 0 || rows.size() != count+1) {
    error("EDSRESTRAINTS: NRES must match the number of rows.");
    return;
  }
  std::set<unsigned int> ids;
  const unsigned int states = sim.param().eds.numstates;
  for (size_t row = 1; row < rows.size(); ++row) {
    std::istringstream line(rows[row]);
    long long index = 0;
    int type = -1;
    simulation::Parameter::colvar_bias_spec spec;
    if (!(line >> index >> type >> spec.method >> spec.rah >> spec.tau >> spec.virial >> spec.linear_tail)
        || index <= 0 || index > std::numeric_limits<unsigned int>::max() || type != 0
        || spec.method != 0 || spec.tau != 0 || spec.virial > 1
        || !std::isfinite(spec.linear_tail) || spec.linear_tail < 0) {
      error("EDSRESTRAINTS: expected positive INDX, TYPE 0 DISTANCE, METHOD 0, RAH, TAU 0, VIRIAL 0/1, finite CVLIN >= 0.");
      continue;
    }
    spec.index = index;
    spec.type = "DISTANCE";
    validate_rah(spec.rah,spec.type,"EDSRESTRAINTS");
    std::vector<simulation::Parameter::colvar_bias_spec> parameters;
    for (unsigned int state = 0; state < states; ++state) {
      if (!(line >> spec.target >> spec.k) || !std::isfinite(spec.target)
          || !std::isfinite(spec.k) || spec.target < 0 || spec.k < 0) break;
      parameters.push_back(spec);
    }
    if (parameters.size() != states || (line >> extra)) {
      error("EDSRESTRAINTS: supply exactly one finite nonnegative TARGET CVK pair per EDS state.");
      continue;
    }
    if (!ids.insert(spec.index).second) {
      error("EDSRESTRAINTS: duplicate INDX.");
      continue;
    }
    sim.param().colvarres.eds_bias_specs.push_back(parameters);
  }
  double dish = 0, disc = 0;
  std::istringstream lengths(geometry.front());
  if (!(lengths >> dish >> disc) || (lengths >> extra)
      || !std::isfinite(dish) || !std::isfinite(disc) || dish <= 0 || disc <= 0) {
    error("EDSDISTANCERESSPEC: expected positive finite DISH DISC.");
    return;
  }
  std::map<unsigned int,std::pair<util::Virtual_Atom,util::Virtual_Atom> > endpoints;
  for (size_t row = 1; row < geometry.size(); ++row) {
    std::istringstream line(geometry[row]);
    long long index = 0;
    if (!(line >> index) || index <= 0 || index > std::numeric_limits<unsigned int>::max()
        || !ids.count(index)) {
      error("EDSDISTANCERESSPEC: INDX must match an EDSRESTRAINTS row.");
      continue;
    }
    std::vector<util::Virtual_Atom> atoms;
    for (int endpoint = 0; endpoint < 2; ++endpoint) {
      std::vector<int> real_atoms;
      bool valid = true;
      for (int slot = 0; slot < 4; ++slot) {
        int atom = 0;
        if (!(line >> atom) || atom < 0 || static_cast<unsigned int>(atom) > topo.num_atoms()) valid = false;
        if (atom > 0 && static_cast<unsigned int>(atom) <= topo.num_atoms()) real_atoms.push_back(atom-1);
      }
      int type = 0;
      if (!(line >> type)) valid = false;
      size_t expected = 0;
      switch (type) {
        case 0: expected = 1; break;
        case 1: expected = 4; break;
        case 2: case 3: case 4: case 6: case 8: expected = 3; break;
        case 5: case 7: expected = 2; break;
        case -1: case -2: expected = real_atoms.size(); break;
        default: valid = false;
      }
      if (!valid || real_atoms.empty() || real_atoms.size() != expected) break;
      atoms.push_back(util::Virtual_Atom(util::virtual_type(type),real_atoms,dish,disc));
    }
    if (atoms.size() != 2 || (line >> extra)) {
      error("EDSDISTANCERESSPEC: expected INDX and two valid i j k l TYPE virtual endpoints.");
      continue;
    }
    if (!endpoints.emplace(index,std::make_pair(atoms[0],atoms[1])).second)
      error("EDSDISTANCERESSPEC: duplicate INDX.");
  }
  for (const auto &parameters : sim.param().colvarres.eds_bias_specs) {
    const auto found = endpoints.find(parameters[0].index);
    if (found == endpoints.end()) {
      error("EDSDISTANCERESSPEC: missing geometry for an EDSRESTRAINTS row.");
      continue;
    }
    std::vector<double> targets, strengths;
    for (const auto &spec : parameters) {
      targets.push_back(spec.target);
      strengths.push_back(spec.k);
    }
    topo.eds_distance_restraints().push_back(topology::eds_distance_restraint_struct(
      found->second.first,found->second.second,targets,strengths,parameters[0].rah));
  }
}

} // namespace colvar_input
} // namespace io
