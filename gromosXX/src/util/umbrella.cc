/**
 * @file umbrella.cc
 * umbrella potentials for LEUS
 */
#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../math/periodicity.h"

#include "../util/le_coordinate.h"
#include "../util/umbrella_weight.h"
#include "../util/template_split.h"
#include "../util/debug.h"

#include "umbrella.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE leus

util::Umbrella::Umbrella(int id, unsigned int dim, Umbrella_Weight_Factory * factory) {
  this->id = id;
  // don't resize the coordinates. They'll get attached
  functional_form.resize(dim);
  variable_type.resize(dim);
  width.resize(dim);
  width_rel.resize(dim);
  cutoff.resize(dim);
  cutoff_rel.resize(dim);
  num_grid_points.resize(dim);
  grid_min.resize(dim);
  grid_min_rel.resize(dim);
  grid_max.resize(dim);
  grid_max_rel.resize(dim);
  grid_spacing_rel.resize(dim);
  building = false;
  enabled = false;
  umbrella_weight_factory = factory;
  configuration_block = "";
  configuration_block_pos = 0;

}

util::Umbrella::Umbrella(const Umbrella& u) {
  *this = u;
}

util::Umbrella & util::Umbrella::operator=(const util::Umbrella& u) {
  id = u.id;
  coordinates = u.coordinates;
  force_constant = u.force_constant;
  functional_form = u.functional_form;
  variable_type = u.variable_type;
  width = u.width;
  width_rel = u.width_rel;
  cutoff = u.cutoff;
  cutoff_rel = u.cutoff_rel;
  num_grid_points = u.num_grid_points;
  grid_min = u.grid_min;
  grid_min_rel = u.grid_min_rel;
  grid_max = u.grid_max;
  grid_max_rel = u.grid_max_rel;
  grid_spacing_rel = u.grid_spacing_rel;
  configurations = u.configurations;
  building = u.building;
  enabled = u.enabled;
  umbrella_weight_factory = u.umbrella_weight_factory;
  configuration_block = u.configuration_block;
  configuration_block_pos = u.configuration_block_pos;

  return *this;
}

util::Umbrella::~Umbrella() {
/*
  if (umbrella_weight_factory != NULL)
    delete umbrella_weight_factory;
*/
}


bool util::Umbrella::leus_conf::operator <(const leus_conf& c) const {
  assert(pos.size() == c.pos.size());
  for(unsigned int i = 0; i < pos.size(); ++i) {
    if (pos[i] < c.pos[i])
      return true;
    else if (pos[i] > c.pos[i])
      return false;
    else
      continue;
  }
  return false;
}

void util::Umbrella::build(
        configuration::Configuration & conf) {
  DEBUG(5, "Building umbrella " << id);
  const unsigned int dim = this->dim();
  const unsigned int num_crd = coordinates.size() / dim;
  for (unsigned int crd = 0; crd < num_crd; ++crd) {
    bool outside = false;
    std::vector<int> pos(dim);
    for (unsigned int i = 0; i < dim; ++i) {
      DEBUG(6, "\tLE coordinate: " << i);

      // get the value and grid it
      const double qi = coordinates[crd * dim + i]->get_value(grid_min_rel[i], grid_max_rel[i]);
      DEBUG(6, "\tqi: " << qi);

      pos[i] = int(floor((qi - grid_min_rel[i]) * grid_spacing_rel[i] + 0.5));
      DEBUG(6, "\tpos: " << pos[i]);
      if (pos[i] < 0 || pos[i] >= int(num_grid_points[i])) {
        DEBUG(8, "\t\toutside the grid");
        outside = true;
        break;
      }
    }

    // don't add it if it is not on the grid!
    if (outside)
      return;

    const util::Umbrella::leus_conf grid_point(pos);
    // binary search this grid point
    std::map<util::Umbrella::leus_conf, util::Umbrella_Weight*>::iterator conf_it =
            configurations.find(grid_point);

    if (conf_it != configurations.end()) {
      // found
      conf_it->second->increment_weight();
    } else {
      // not found so add the grid point to the visited points
      configurations.insert(
              std::pair<util::Umbrella::leus_conf, util::Umbrella_Weight* > (
              grid_point, umbrella_weight_factory->get_instance()));
    }
  } // coordinate sets
}

void util::Umbrella::apply(
        configuration::Configuration & conf) {
  DEBUG(6, "Applying umbrella " << id);
  const unsigned int dim = this->dim();

  // loop over the coordinate sets
  const unsigned int num_crd = coordinates.size() / dim;
  DEBUG(6, "dimensionality of umbrella: " << dim);
  DEBUG(6, "number of coordinate sets: " << num_crd);
  DEBUG(6, "number of coordinates: " << coordinates.size());
  for (unsigned int crd = 0; crd < num_crd; ++crd) {
    // loop over the stored configurations
    std::map<util::Umbrella::leus_conf, util::Umbrella_Weight*>::const_iterator conf_it =
            configurations.begin(), conf_to = configurations.end();

    for (; conf_it != conf_to; ++conf_it) {
      bool outside = false;
      double energy = 1.0;
      std::vector<double> deriv(dim, 1.0);
      // loop over the coordinates attached to the umbrella
      for (unsigned int i = 0; i < dim; ++i) {
        DEBUG(8, "\tpos[" << i + 1 << "]: " << conf_it->first.pos[i]);
        // get the value of the grid point
        const double qi_0 = grid_min_rel[i] + conf_it->first.pos[i] / grid_spacing_rel[i];
        // get the deviation from the value
        const double dq = coordinates[crd * dim + i]->get_deviation(qi_0);
        DEBUG(8, "\tqi_0: " << qi_0 << " deltaq: " << dq);

        // check cutoff
        const double dq_abs = fabs(dq);
        if (dq_abs > cutoff_rel[i]) {
          DEBUG(8, "\t\tout of cutoff.");
          outside = true;
          break;
        }
        DEBUG(8, "\t\twithin cutoff.");

        // calculate g function and the derivative by Qi
        double g = 0.0, dgdQ = 0.0;
        switch (functional_form[i]) {
          case util::Umbrella::ff_polynomial:
          {
            const double term1 = dq_abs / (width_rel[i] * width_rel[i]);
            const double term2 = term1 * dq_abs / width_rel[i];
            g = 1.0 - 3.0 * term1 * dq_abs + 2.0 * term2 * dq_abs;
            const double sign = dq < 0.0 ? -1.0 : 1.0;
            dgdQ = sign * 6.0 * (-term1 + term2);
            break;
          }
          case util::Umbrella::ff_gaussian:
          {
            const double factor = -1.0 / (width_rel[i] * width_rel[i]);
            const double exp_term = exp(0.5 * dq * dq * factor);
            g = exp_term;
            dgdQ = factor * exp_term * dq;
            break;
          }
          default:
            io::messages.add("Functional form not implemented",
                    "Umbrella",
                    io::message::critical);
        }

        DEBUG(5, "\tg function: " << g);
        DEBUG(5, "\tdgdQ      : " << dgdQ);
        // multiply the g functions
        energy *= g;

        // this is to take care of the product rule
        // derivative is the product of all g functions except the one i=j
        // where it is the derivative of g
        for (unsigned int j = 0; j < dim; ++j) {
          if (j == i)
            deriv[j] *= dgdQ;
          else
            deriv[j] *= g;
        }
      }

      if (outside) {
        DEBUG(8, "\tconf is out of cutoff.");
        continue;
      }

      energy *= force_constant * conf_it->second->get_weight();
      DEBUG(4, "\tenergy: " << energy)
      // store the energy, apply the force
      conf.current().energies.leus_total += energy;
      for (unsigned int i = 0; i < dim; ++i) {
        const double d = deriv[i] * force_constant * conf_it->second->get_weight();
        DEBUG(7, "\tderiv[" << i << "]: " << d);
        coordinates[crd * dim + i]->apply(d);
      }
    } // for visited configurations
  } // for coordinate sets
  DEBUG(6,"\tLEUS energy: " << conf.current().energies.leus_total);
}

void util::Umbrella::transform_units() {
  for(unsigned int i = 0; i < dim(); ++i) {
    double conversion_factor = 1.0;
    double conversion_factor_grid = 1.0;
    switch(variable_type[i]) {
      case util::Umbrella::vt_dihedral:
        conversion_factor = math::Pi / 180.0;
        conversion_factor_grid = math::Pi * 2.0 / num_grid_points[i];
        break;
      default:
        conversion_factor_grid = (grid_max[i] - grid_min[i])/(num_grid_points[i]-1.0);
        break;
    }
    cutoff_rel[i] = conversion_factor_grid * cutoff[i];
    width_rel[i] = conversion_factor_grid * width[i];
    grid_min_rel[i] = conversion_factor * grid_min[i];
    grid_max_rel[i] = conversion_factor * grid_max[i];

    // check for weird periodicity settings
    if (grid_max_rel[i] != grid_min_rel[i])
      grid_spacing_rel[i] = (num_grid_points[i]-1.0) / (grid_max_rel[i] - grid_min_rel[i]);
    else {
      switch(variable_type[i]) {
        case util::Umbrella::vt_dihedral :
          grid_spacing_rel[i] = num_grid_points[i] / (2.0 * math::Pi);
          grid_max_rel[i] = grid_min_rel[i] + 2.0 * math::Pi;
          break;
        default:
          io::messages.add("variable type cannot be periodic", "Umbrella", io::message::error);
      }
    }
  }
}

void util::Umbrella::calculate_coordinates(
        configuration::Configuration & conf) {
  DEBUG(4, "Calculating coordinates of umbrella " << id);
  for (unsigned int i = 0; i < coordinates.size(); ++i) {
    coordinates[i]->calculate(conf);
  }
}

void util::Umbrella::read_configuration() {
  DEBUG(4, "Converting read configurations of umbrella " << id);
  std::istringstream _lineStream(configuration_block);
  _lineStream.seekg(configuration_block_pos);

  unsigned int nconle = 0;
  _lineStream >> nconle;
  if (_lineStream.fail()) {
    io::messages.add("LEUSBIAS block: Could not read umbrella number of configurations",
            "Umbrella", io::message::error);
    return;
  }
  // loop over sampled points
  for (unsigned int n = 0; n < nconle; ++n) {
    util::Umbrella::leus_conf cnf(dim());
    util::Umbrella_Weight * weight = umbrella_weight_factory->get_instance();
    _lineStream >> (*weight);
    for (unsigned int i = 0; i < dim(); ++i) {
      _lineStream >> cnf.pos[i];
      if (_lineStream.fail()) {
        io::messages.add("LEUSBIAS block: Could not read stored configurations",
                "Umbrella", io::message::error);
        return;
      }
      --cnf.pos[i]; // our arrays start at 0 and not 1 as in the format
    }
    configurations.insert(std::pair<util::Umbrella::leus_conf, util::Umbrella_Weight*>(cnf, weight));
  } // for configurations
}

std::string util::Umbrella::str() const {
  std::ostringstream os;
  os << "  Umbrella (" << id << "): ";
  if (enabled) {
    os << "enabled, " << (building ? "building" : "frozen") << ".";
  } else {
    os << "disabled.";
  }
  os.precision(8);
  os << "\n    force constant: " << std::setw(15) << force_constant;
  os << "\n    " << configurations.size() << " configurations visited.\n";
  os << "\n    dimensionality: " << dim();
  os << "\n    functional form/grid properties:\n";
  for (unsigned int i = 0; i < dim(); ++i) {
    os << "    - function    : ";
    switch (functional_form[i]) {
      case ff_gaussian:
        os << "gaussian";
        break;
      case ff_polynomial:
        os << "truncated polynomial";
        break;
      default:
        os << "unkown";
    }
    os << "\n";
    os << "      width       : " << std::setw(15) << width[i] << "\n"
       << "      cutoff      : " << std::setw(15) << cutoff[i] << "\n"
       << "      grid points : " << num_grid_points[i] << "\n"
       << "      grid min    : " << std::setw(15) << grid_min[i] << "\n"
       << "      grid max    : " << std::setw(15) << grid_max[i] << "\n";
  }
  unsigned int num_crd = coordinates.size() / dim();
  os << "    " << num_crd << " coordinate sets:\n";
  for (unsigned int crd = 0; crd < num_crd; ++crd) {
    os << "    Coordinate set: " << crd+1 << "\n";
    for (unsigned int i = 0; i < dim(); ++i) {
      os << "    - " << coordinates[crd*dim()+i]->str() << "\n";
    }
  }

  return os.str();
}

