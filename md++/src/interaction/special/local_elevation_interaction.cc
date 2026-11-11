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
 * @file local_elevation_interaction.cc
 * apply LEUS
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/special/local_elevation_interaction.h"
#include "../../util/le_coordinate.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"
#include <vector>
#include <map>

#include "util/umbrella_weight.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

int interaction::Local_Elevation_Interaction::
calculate_interactions(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim) {
  DEBUG(4, "local elevation interaction");
  m_timer.start(sim);
  auto &umbrellas = conf.special().umbrellas;
  std::vector<int> deposit(umbrellas.size(),0);
  // Validate every geometry and guard before changing any history.
  for (size_t i=0;i<umbrellas.size();++i) {
    auto &umbrella=umbrellas[i];
    if (!umbrella.enabled) continue;
    if (umbrella.calculate_coordinates(topo,conf,sim)) { m_timer.stop(); return 1; }
    if (umbrella.building) {
      deposit[i]=umbrella.deposition_due(sim.time(),sim.time_step_size());
      if (deposit[i]<0) { m_timer.stop(); return 1; }
    }
  }
  for (size_t i=0;i<umbrellas.size();++i) {
    auto &umbrella=umbrellas[i];
    if (!umbrella.enabled) continue;
    if (deposit[i]) {
      umbrella.build(conf);
      umbrella.mark_deposition(sim.time(),sim.time_step_size());
    }
    umbrella.apply(conf,sim.param().localelev.virial!=0);
    if (io::messages.contains(io::message::error) || io::messages.contains(io::message::critical)) {
      m_timer.stop(); return 1;
    }
  }

  m_timer.stop();
  return 0;
}

int interaction::Local_Elevation_Interaction::init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os,
            bool quiet) {
  std::set<int> umbrella_ids;
  for (const auto &umbrella : conf.special().umbrellas)
    if (!umbrella_ids.insert(umbrella.id).second) {
      io::messages.add("Duplicate LE umbrella ID.","Local_Elevation_Interaction",io::message::error);
      return 1;
    }
  // Rebind to this topology, including after cloning/reinitialisation.
  for (auto &umbrella : conf.special().umbrellas) {
    umbrella.coordinates.clear();
    umbrella.enabled=false;
    umbrella.building=false;
  }
  // attach coordinates to umbrellas
  std::vector<util::LE_Coordinate*>::const_iterator le_it = topo.le_coordinates().begin(),
          le_to = topo.le_coordinates().end();
  for(; le_it != le_to; ++le_it) {
    // find the umbrella
    bool found = false;
    std::vector<util::Umbrella>::iterator umb_it = conf.special().umbrellas.begin(),
            umb_to = conf.special().umbrellas.end();
    for(; umb_it != umb_to; ++umb_it) {
      if (umb_it->id == (*le_it)->umbrella_id()) {
        umb_it->coordinates.push_back(*le_it);
        found = true; break;
      }
    }
    if (!found) {
      std::ostringstream msg;
      msg << "LE coordinate points to unkown umbrella id ("
              << (*le_it)->umbrella_id() << ")";
      io::messages.add(msg.str(), "Local_Elevation_Interaction", io::message::error);
      return 1;
    }
  }

  // enable the umbrella if required. decide whether we build or freeze
  std::map<int, bool>::const_iterator id_it = sim.param().localelev.umbrellas.begin(),
          id_to = sim.param().localelev.umbrellas.end();
  for(; id_it != id_to; ++id_it) {
    // find the umbrella
    bool found = false;
    std::vector<util::Umbrella>::iterator umb_it = conf.special().umbrellas.begin(),
            umb_to = conf.special().umbrellas.end();
    for(; umb_it != umb_to; ++umb_it) {
      if (umb_it->id == id_it->first) {
        umb_it->enabled = true;
        umb_it->building = id_it->second;
        found = true; break;
      }
    }
    if (!found) {
      std::ostringstream msg;
      msg << "Could not find umbrella (" << id_it->first << ") in database";
      io::messages.add(msg.str(), "Local_Elevation_Interaction", io::message::error);
      return 1;
    }
  }

  // check for problems
  std::vector<util::Umbrella>::iterator umb_it = conf.special().umbrellas.begin(),
          umb_to = conf.special().umbrellas.end();
  for (; umb_it != umb_to; ++umb_it) {
    // attach weight to umbrella if this was not done yet
    umb_it->initialise_weight_factory();
    if (!umb_it->validate_grid()) return 1;
    if (umb_it->enabled && sim.param().localelev.virial) {
      if (!dynamic_cast<util::Number_Of_Visits_Umbrella_Weight_Factory*>(umb_it->umbrella_weight_factory)) {
        io::messages.add("LE virial is currently supported only with visit-count weights.",
                         "Local_Elevation_Interaction",io::message::error);
        return 1;
      }
      for (auto coordinate : umb_it->coordinates) {
        if (!coordinate->supports_virial()) {
          io::messages.add("LE virial requires CV-backed coordinates; distance-field/legacy coordinates are unsupported.",
                           "Local_Elevation_Interaction",io::message::error);
          return 1;
        }
      }
    }

    // dimensions
    if (umb_it->enabled && umb_it->coordinates.empty()) {
      io::messages.add("Enabled LE umbrella has no attached coordinates.","Local_Elevation_Interaction",io::message::error);
      return 1;
    }
    const unsigned int num_crd = umb_it->coordinates.size() / umb_it->dim();
    if (umb_it->coordinates.size() % umb_it->dim() != 0) {
      std::ostringstream msg;
      msg << "Umbrella (" << umb_it->id << ") has dimension "
          << umb_it->variable_type.size() << " but " << umb_it->dim()
          << " variables attached to it. The number of variables attached has to be a multiple of the dimension.";
      io::messages.add(msg.str(), "Local_Elevation_Interaction", io::message::error);
      return 1;
    }

    // types
    for (unsigned int crd = 0; crd < num_crd; ++crd) {
      for (unsigned int i = 0; i < umb_it->dim(); ++i) {
        int vt = umb_it->variable_type[i];
        unsigned int crd_index = i + crd*umb_it->dim();
        int le_vt = umb_it->coordinates[crd_index]->get_type();
        if (vt != le_vt) {
          std::ostringstream msg;
          msg << "Umbrella (" << umb_it->id << ") LE coordinate " << (crd_index + 1)
                  << " variable types do not match: " << vt << " vs. " << le_vt << ".";
          io::messages.add(msg.str(), "Local_Elevation_Interaction", io::message::error);
        }
      }
    }

    const std::string signature=umb_it->definition_signature();
    if (!umb_it->restart_signature.empty() && umb_it->restart_signature!=signature) {
      io::messages.add("LEUSSTATE coordinate/grid mismatch for umbrella "+std::to_string(umb_it->id),
                      "Local_Elevation_Interaction",io::message::error);
      return 1;
    }
    umb_it->restart_signature=signature;
    if (umb_it->building && umb_it->has_deposition
        && umb_it->deposition_due(sim.time(),sim.time_step_size())<0) return 1;
    // transform the units
    umb_it->transform_units();
    // read the configurations
    if (umb_it->configurations.empty()) umb_it->read_configuration();
    if (io::messages.contains(io::message::error) || io::messages.contains(io::message::critical)) return 1;
  }

  if (!quiet) {
    os << "LE bias virial: " << (sim.param().localelev.virial ? "on" : "off") << '\n';
    os << "LOCALELEV\n"
            << "Local-elevation umbrella sampling is on.\n" 
            << "Number of umbrellas: " << conf.special().umbrellas.size() << "\n\n";
    std::vector<util::Umbrella>::const_iterator umb_it = conf.special().umbrellas.begin(),
            umb_to = conf.special().umbrellas.end();
    for(; umb_it != umb_to; ++umb_it) {
      os << umb_it->str() << "\n";
    }
    os << "END\n";
  }

  return 0;
}
