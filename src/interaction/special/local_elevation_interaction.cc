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
  m_timer.start();
  // loop over umbrellas
  std::vector<util::Umbrella>::iterator
  it = conf.special().umbrellas.begin(), to = conf.special().umbrellas.end();
  for(; it != to; ++it) {
    if (it->enabled) {
      it->calculate_coordinates(conf);
      // built and/or apply it
      if (it->building) {
        it->build(conf);
      }
      it->apply(conf);
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
    if (umb_it->umbrella_weight_factory == NULL)
      umb_it->umbrella_weight_factory = new util::Number_Of_Visits_Umbrella_Weight_Factory();

    // dimensions
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

    // transform the units
    umb_it->transform_units();
    // read the configurations
    umb_it->read_configuration();
  }

  if (!quiet) {
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
