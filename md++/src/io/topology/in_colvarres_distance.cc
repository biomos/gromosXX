/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres_distance.cc
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

void io::In_Colvarres::read_distanceresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

  // -------------------------------------------------------------------------
  // DISTANCERESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "DISTANCERESSPEC";
    DEBUG(10, "DISTANCERESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim, "DISTANCE") > 0) {
        io::messages.add("RESTRAINTS contains DISTANCE entries but no non-empty "
                         "DISTANCERESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);

      std::vector<std::string>::const_iterator it = buffer.begin() + 1;
      const std::vector<std::string>::const_iterator to = buffer.end() - 1;

      for (; it != to && !line_is_data(*it); ++it) {}

      double dish = 0.0, disc = 0.0;
      if (it != to) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> dish >> disc;
      }
      if (_lineStream.fail()) {
        io::messages.add("bad first data line in DISTANCERESSPEC block: expected DISH DISC",
                         "In_Colvarres", io::message::error);
      }

      if (!quiet) {
        os << std::setw(10) << "DISH"
           << std::setw(10) << "DISC"
           << "\n"
           << std::setw(10) << dish
           << std::setw(10) << disc
           << "\n";
      }

      size_t distance_index = 0;
      bool indexed_distance_geometries = false;
      std::vector<std::pair<unsigned int, topology::distance_restraint_struct> >
        indexed_distances;
      if (it != to) ++it;
      for (; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        int type1 = 0, type2 = 0;
        unsigned int geom_index = 0;
        std::vector<int> atom1, atom2;
        std::vector<std::string> tokens = tokens_from_line(*it);

        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 11) {
          _lineStream >> geom_index;
        }

        for (unsigned int i = 0; i < io::In_Colvarres::MAX_ATOMS; ++i) {
          unsigned int atom = 0;
          _lineStream >> atom;
          if (atom > topo.num_atoms()) {
            std::ostringstream msg;
            msg << blockname << " block: atom number out of range: "
                << atom << ", last atom is " << topo.num_atoms();
            io::messages.add(msg.str(), "In_Colvarres", io::message::error);
          }
          if (atom > 0) atom1.push_back(atom - 1);
        }
        _lineStream >> type1;

        for (unsigned int i = 0; i < io::In_Colvarres::MAX_ATOMS; ++i) {
          unsigned int atom = 0;
          _lineStream >> atom;
          if (atom > topo.num_atoms()) {
            std::ostringstream msg;
            msg << blockname << " block: atom number out of range: "
                << atom << ", last atom is " << topo.num_atoms();
            io::messages.add(msg.str(), "In_Colvarres", io::message::error);
          }
          if (atom > 0) atom2.push_back(atom - 1);
        }
        _lineStream >> type2;

        if (_lineStream.fail()) {
          io::messages.add("bad line in DISTANCERESSPEC block: expected two virtual atom definitions",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const simulation::Parameter::colvar_bias_spec *spec =
          geom_index > 0
            ? indexed_colvar_spec(sim.param().colvarres.bias_specs,
                                  "DISTANCE", geom_index)
            : nth_colvar_spec(sim, "DISTANCE", distance_index);
        if (spec == NULL) {
          if (geom_index > 0) {
            std::ostringstream msg;
            msg << "DISTANCERESSPEC index " << geom_index
                << " has no matching RESTRAINTS DISTANCE entry; geometry is ignored.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
          else {
            io::messages.add("DISTANCERESSPEC contains more geometries than RESTRAINTS DISTANCE entries; extra geometries are ignored.",
                             "In_Colvarres", io::message::warning);
          }
          ++distance_index;
          continue;
        }

        util::Virtual_Atom v1(util::virtual_type(type1), atom1, dish, disc);
        util::Virtual_Atom v2(util::virtual_type(type2), atom2, dish, disc);

        topology::distance_restraint_struct restraint(v1, v2,
                                                      spec->target, 1.0,
                                                      spec->rah);
        if (geom_index > 0) {
          indexed_distance_geometries = true;
          indexed_distances.push_back(std::make_pair(geom_index, restraint));
        }
        else {
          topo.distance_restraints().push_back(restraint);
        }

        ++distance_index;
      }

      if (indexed_distance_geometries) {
        for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
             spec_it = sim.param().colvarres.bias_specs.begin(),
             spec_to = sim.param().colvarres.bias_specs.end();
             spec_it != spec_to; ++spec_it) {
          if (spec_it->type != "DISTANCE" || spec_it->index == 0) continue;
          bool found = false;
          for (std::vector<std::pair<unsigned int, topology::distance_restraint_struct> >::const_iterator
               geom_it = indexed_distances.begin(),
               geom_to = indexed_distances.end(); geom_it != geom_to; ++geom_it) {
            if (geom_it->first != spec_it->index) continue;
            topo.distance_restraints().push_back(geom_it->second);
            found = true;
            break;
          }
          if (!found) {
            std::ostringstream msg;
            msg << "RESTRAINTS DISTANCE index " << spec_it->index
                << " has no matching DISTANCERESSPEC geometry.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
        }
      }
    }
  }
}

void io::In_Colvarres::read_pertdistanceresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

  // -------------------------------------------------------------------------
  // PERTDISTANCERESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "PERTDISTANCERESSPEC";
    DEBUG(10, "PERTDISTANCERESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim.param().colvarres.pert_bias_specs, "DISTANCE") > 0) {
        io::messages.add("PERTRESTRAINTS contains DISTANCE entries but no non-empty "
                         "PERTDISTANCERESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);

      std::vector<std::string>::const_iterator it = buffer.begin() + 1;
      const std::vector<std::string>::const_iterator to = buffer.end() - 1;
      for (; it != to && !line_is_data(*it); ++it) {}

      double dish = 0.0, disc = 0.0;
      if (it != to) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> dish >> disc;
      }
      if (_lineStream.fail()) {
        io::messages.add("bad first data line in PERTDISTANCERESSPEC block: expected DISH DISC",
                         "In_Colvarres", io::message::error);
      }

      size_t distance_index = 0;
      bool indexed_distance_geometries = false;
      std::vector<std::pair<unsigned int, topology::distance_restraint_struct> >
        indexed_distances;
      if (it != to) ++it;
      for (; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        int type1 = 0, type2 = 0;
        unsigned int geom_index = 0;
        std::vector<int> atom1, atom2;
        std::vector<std::string> tokens = tokens_from_line(*it);
        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 11) {
          _lineStream >> geom_index;
        }

        for (unsigned int i = 0; i < io::In_Colvarres::MAX_ATOMS; ++i) {
          unsigned int atom = 0;
          _lineStream >> atom;
          if (atom > topo.num_atoms()) {
            io::messages.add("PERTDISTANCERESSPEC block: atom number out of range",
                             "In_Colvarres", io::message::error);
          }
          if (atom > 0) atom1.push_back(atom - 1);
        }
        _lineStream >> type1;

        for (unsigned int i = 0; i < io::In_Colvarres::MAX_ATOMS; ++i) {
          unsigned int atom = 0;
          _lineStream >> atom;
          if (atom > topo.num_atoms()) {
            io::messages.add("PERTDISTANCERESSPEC block: atom number out of range",
                             "In_Colvarres", io::message::error);
          }
          if (atom > 0) atom2.push_back(atom - 1);
        }
        _lineStream >> type2;

        if (_lineStream.fail()) {
          io::messages.add("bad line in PERTDISTANCERESSPEC block: expected two virtual atom definitions",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const simulation::Parameter::colvar_bias_spec *spec =
          geom_index > 0
            ? indexed_colvar_spec(sim.param().colvarres.pert_bias_specs,
                                  "DISTANCE", geom_index)
            : nth_colvar_spec(sim.param().colvarres.pert_bias_specs,
                              "DISTANCE", 2 * distance_index);
        if (spec == NULL) {
          if (geom_index > 0) {
            std::ostringstream msg;
            msg << "PERTDISTANCERESSPEC index " << geom_index
                << " has no matching PERTRESTRAINTS DISTANCE entry; geometry is ignored.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
          else {
            io::messages.add("PERTDISTANCERESSPEC contains more geometries than PERTRESTRAINTS DISTANCE entries; extra geometries are ignored.",
                             "In_Colvarres", io::message::warning);
          }
          ++distance_index;
          continue;
        }

        util::Virtual_Atom v1(util::virtual_type(type1), atom1, dish, disc);
        util::Virtual_Atom v2(util::virtual_type(type2), atom2, dish, disc);
        topology::distance_restraint_struct restraint(v1, v2,
                                                      spec->target, 1.0,
                                                      spec->rah);
        if (geom_index > 0) {
          indexed_distance_geometries = true;
          indexed_distances.push_back(std::make_pair(geom_index, restraint));
        }
        else {
          topo.distance_restraints().push_back(restraint);
        }
        ++distance_index;
      }

      if (indexed_distance_geometries) {
        std::set<unsigned int> appended;
        for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
             spec_it = sim.param().colvarres.pert_bias_specs.begin(),
             spec_to = sim.param().colvarres.pert_bias_specs.end();
             spec_it != spec_to; ++spec_it) {
          if (spec_it->type != "DISTANCE" || spec_it->index == 0
              || appended.count(spec_it->index)) continue;
          bool found = false;
          for (std::vector<std::pair<unsigned int, topology::distance_restraint_struct> >::const_iterator
               geom_it = indexed_distances.begin(),
               geom_to = indexed_distances.end(); geom_it != geom_to; ++geom_it) {
            if (geom_it->first != spec_it->index) continue;
            topo.distance_restraints().push_back(geom_it->second);
            appended.insert(spec_it->index);
            found = true;
            break;
          }
          if (!found) {
            std::ostringstream msg;
            msg << "PERTRESTRAINTS DISTANCE index " << spec_it->index
                << " has no matching PERTDISTANCERESSPEC geometry.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
        }
      }
    }
  }
}

namespace io {
namespace colvar_input {

void read_distance_difference_geometries(const std::vector<std::string> &buffer,
    const std::string &blockname,
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    bool perturbed, topology::Topology &topo, std::set<std::string> &block_read)
{
  const std::string type = "DISTANCE_DIFFERENCE";
  if (buffer.size() <= 2) {
    if (count_colvar_specs(specs, type) != 0)
      io::messages.add("Missing " + blockname + " for DISTANCE_DIFFERENCE biases.",
                       "In_Colvarres", io::message::error);
    return;
  }
  block_read.insert(blockname);
  std::map<unsigned int, topology::distance_difference_restraint_struct> geometries;
  double dish = 0.0, disc = 0.0;
  bool lengths_read = false;
  for (auto it = buffer.begin() + 1; it != buffer.end() - 1; ++it) {
    if (!line_is_data(*it)) continue;
    std::istringstream line(it->substr(0, it->find('#')));
    std::string extra;
    if (!lengths_read) {
      lengths_read = true;
      if (!(line >> dish >> disc) || (line >> extra)
          || !std::isfinite(dish) || !std::isfinite(disc) || dish <= 0 || disc <= 0) {
        io::messages.add(blockname + ": expected positive finite DISH DISC.",
                         "In_Colvarres", io::message::error);
        return;
      }
      continue;
    }
    long long input_index = 0;
    if (!(line >> input_index) || input_index <= 0
        || input_index > std::numeric_limits<unsigned int>::max()) {
      io::messages.add(blockname + ": expected a positive restraint index.",
                       "In_Colvarres", io::message::error);
      continue;
    }
    const unsigned int index = static_cast<unsigned int>(input_index);
    std::vector<util::Virtual_Atom> endpoints;
    for (int endpoint = 0; endpoint < 4; ++endpoint) {
      std::vector<int> atoms;
      bool valid = true;
      for (int slot = 0; slot < 4; ++slot) {
        unsigned int atom = 0;
        if (!(line >> atom) || atom > topo.num_atoms()) valid = false;
        if (atom > 0 && atom <= topo.num_atoms()) atoms.push_back(atom - 1);
      }
      int va_type = 0;
      if (!(line >> va_type)) valid = false;
      size_t expected = 0;
      switch (va_type) {
        case 0: expected = 1; break;
        case 1: expected = 4; break;
        case 2: case 3: case 4: case 6: case 8: expected = 3; break;
        case 5: case 7: expected = 2; break;
        case -1: case -2: expected = atoms.size(); break;
        default: valid = false;
      }
      if (!valid || atoms.empty() || atoms.size() != expected) break;
      endpoints.push_back(util::Virtual_Atom(util::virtual_type(va_type), atoms, dish, disc));
    }
    if (endpoints.size() != 4 || (line >> extra)) {
      io::messages.add(blockname + ": expected INDX followed by four valid i j k l TYPE endpoints; check atom numbers and site sizes.",
                       "In_Colvarres", io::message::error);
      continue;
    }
    if (indexed_colvar_spec(specs, type, index) == NULL) {
      io::messages.add(blockname + ": geometry index has no matching DISTANCE_DIFFERENCE bias.",
                       "In_Colvarres", io::message::error);
      continue;
    }
    if (!geometries.insert(std::make_pair(index,
        topology::distance_difference_restraint_struct(endpoints[0], endpoints[1],
                                                       endpoints[2], endpoints[3]))).second)
      io::messages.add(blockname + ": duplicate geometry index.",
                       "In_Colvarres", io::message::error);
  }
  // Match input IDs, then store geometry in bias order, just like the other CVs.
  size_t state = 0;
  for (const auto &spec : specs) {
    if (spec.type != type) continue;
    if (perturbed && state++ % 2 != 0) continue;
    const auto geometry = geometries.find(spec.index);
    if (geometry == geometries.end()) {
      io::messages.add(blockname + ": DISTANCE_DIFFERENCE bias has no matching geometry.",
                       "In_Colvarres", io::message::error);
    } else {
      topo.distance_difference_restraints().push_back(geometry->second);
    }
  }
}

} // namespace colvar_input
} // namespace io
