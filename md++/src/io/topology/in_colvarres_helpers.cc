/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres_helpers.cc
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

bool line_is_data(const std::string &line) {
  const std::string::size_type pos = line.find_first_not_of(" \t\r\n");
  return pos != std::string::npos && line[pos] != '#';
}

std::string colvar_type_from_int(int type) {
  switch (type) {
    case 0: return "DISTANCE";
    case 1: return "ANGLE";
    case 2: return "DIHEDRAL";
    case 3: return "COORDNUM";
    case 4: return "DISTANCE_DIFFERENCE";
    default: return "";
  }
}

void validate_rah(int rah, const std::string &type,
                  const std::string &blockname) {
  // Wall behavior applies to every scalar CV. Projection codes additionally
  // change distance geometry and have no meaning for the other CV types.
  if (rah >= -1 && rah <= 1) return;

  if (type == "DISTANCE") {
    for (int projection = 10; projection <= 60; projection += 10) {
      if (rah >= projection - 1 && rah <= projection + 1) return;
    }
  }

  std::ostringstream msg;
  msg << blockname << " block: invalid RAH " << rah << " for " << type
      << ". Use -1 (lower wall), 0 (full harmonic), or 1 (upper wall)";
  if (type == "DISTANCE") {
    msg << ", optionally added to a distance projection code "
        << "10, 20, 30, 40, 50, or 60";
  }
  msg << ".";
  io::messages.add(msg.str(), "In_Colvarres", io::message::error);
}

std::vector<std::string> tokens_from_line(const std::string &line) {
  std::vector<std::string> tokens;
  std::istringstream stream(line);
  std::string token;
  while (stream >> token) tokens.push_back(token);
  return tokens;
}

bool decode_method(int raw,
                   int &averaging,
                   int &force_scale,
                   bool &constraining,
                   const std::string &type,
                   const std::string &blockname) {
  constraining = false;

  if (raw == -1) {
    averaging = 0;
    force_scale = 0;
    constraining = true;
    io::messages.add(blockname + " block: METHOD -1 colvar constraints are "
                     "recognised but not implemented yet.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  averaging = raw % 10;
  force_scale = raw / 10;

  if (raw < 0 || averaging < 0 || averaging > 2
      || force_scale < 0 || force_scale > 2) {
    io::messages.add(blockname + " block: METHOD must be one of -1, 0, 1, 2, 11, 12, or 22.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale > 0 && averaging == 0) {
    io::messages.add(blockname + " block: force scaling in METHOD requires time averaging.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale == 2 && averaging != 2) {
    io::messages.add(blockname + " block: chain-rule force scaling is only meaningful for METHOD 22.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (type == "DISTANCE_DIFFERENCE" && averaging == 2) {
    io::messages.add(blockname + " block: inverse-cubic averaging is not defined for signed DISTANCE_DIFFERENCE values.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale > 0 && type != "DISTANCE"
      && !(type == "DISTANCE_DIFFERENCE" && raw == 11)) {
    io::messages.add(blockname + " block: force-scaled averaging requires DISTANCE, or METHOD 11 for DISTANCE_DIFFERENCE.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  return true;
}

const simulation::Parameter::colvar_bias_spec *indexed_colvar_spec(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type,
    unsigned int index) {
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->type == type && it->index == index) return &(*it);
  }
  return NULL;
}

bool colvar_index_exists(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    unsigned int index) {
  if (index == 0) return false;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->index == index) return true;
  }
  return false;
}

const simulation::Parameter::colvar_bias_spec *nth_colvar_spec(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type,
    size_t index) {
  size_t n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->type == type) {
      if (n == index) return &(*it);
      ++n;
    }
  }
  return NULL;
}

const simulation::Parameter::colvar_bias_spec *nth_colvar_spec(
    const simulation::Simulation &sim,
    const std::string &type,
    size_t index) {
  return nth_colvar_spec(sim.param().colvarres.bias_specs, type, index);
}

unsigned int count_colvar_specs(const simulation::Simulation &sim,
                                const std::string &type) {
  unsigned int n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = sim.param().colvarres.bias_specs.begin(),
       to = sim.param().colvarres.bias_specs.end(); it != to; ++it) {
    if (it->type == type) ++n;
  }
  return n;
}

unsigned int count_colvar_specs(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type) {
  unsigned int n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->type == type) ++n;
  }
  return n;
}

} // namespace colvar_input
} // namespace io
