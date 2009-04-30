/**
 * @file umbrella.cc
 * umbrella potentials for LEUS
 */
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <math/periodicity.h>

#include <util/le_coordinate.h>
#include <util/template_split.h>
#include <util/debug.h>

#include "umbrella.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE leus

util::Umbrella::Umbrella(int id, unsigned int dim) {
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

  return *this;
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
  DEBUG(1, "Building umbrella " << id);
  const unsigned int dim = coordinates.size();

  bool outside = false;
  std::vector<int> pos(dim);
  for (unsigned int i = 0; i < dim; ++i) {
    DEBUG(1, "\tLE coordinate: " << i);
    
    // get the value and grid it
    const double qi = coordinates[i]->get_value(grid_min_rel[i], grid_max_rel[i]);
    DEBUG(1, "\tqi: " << qi);
    
    pos[i] = floor((qi - grid_min_rel[i]) * grid_spacing_rel[i] + 0.5);
    DEBUG(1, "\tpos: " << pos[i]);
    if (pos[i] < 0 || pos[i] >= int(num_grid_points[i])) {
      DEBUG(1,"\t\toutside the grid");
      outside = true;
      break;
    }
  }

  // don't add it if it is not on the grid!
  if (outside)
    return;

  const util::Umbrella::leus_conf grid_point(pos);
  // binary search this grid point
  std::map<util::Umbrella::leus_conf, unsigned int>::iterator conf_it =
          configurations.find(grid_point);

  if (conf_it != configurations.end()) {
    // found
    conf_it->second++;
  } else {
    // not found so add the grid point to the visited points
    configurations[grid_point] = 1;
  }
}

void util::Umbrella::apply(
        configuration::Configuration & conf) {
  DEBUG(6, "Applying umbrella " << id);

  double & local_elevation_energy = conf.current().energies.leus_total;

  const unsigned int dim = coordinates.size();

  // loop over the stored configurations
  std::map<util::Umbrella::leus_conf, unsigned int>::const_iterator conf_it =
          configurations.begin(), conf_to = configurations.end();

  double energy = 1.0;
  std::vector<double> deriv(dim, 1.0);
  bool outside = false;
  for(; conf_it != conf_to; ++conf_it) {
    // loop over the coordinates attached to the umbrella
    for(unsigned int i = 0; i < dim; ++i) {
      DEBUG(8, "\tpos[" << i+1 << "]: " << conf_it->first.pos[i]);
      // get the value of the grid point
      const double qi_0 = grid_min_rel[i] + conf_it->first.pos[i] / grid_spacing_rel[i];
      // get the deviation from the value
      const double dq = coordinates[i]->get_deviation(qi_0);
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
      double g, dgdQ;
      switch (functional_form[i]) {
        case util::Umbrella::ff_polynomial :
        {
          const double term1 = dq_abs / (width_rel[i] * width_rel[i]);
          const double term2 = term1 * dq_abs / width_rel[i];
          g = 1.0 - 3.0 * term1 * dq_abs + 2.0 * term2 * dq_abs;
          dgdQ = 6.0 * (-term1 + term2);
          break;
        }
        case util::Umbrella::ff_gaussian :
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

      DEBUG(1, "\tg function: " << g);
      // multiply the g functions
      energy *= g;

      // this is to take care of the product rule
      // derivative is the product of all g functions except the one i=j
      // where it is the derivative of g
      for(unsigned int j = 0; j < dim; ++j) {
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

    energy *= force_constant * conf_it->second;
    DEBUG(1,"\tenergy: " << energy)
    // store the energy, apply the force
    local_elevation_energy += energy;
    for(unsigned int i = 0; i < dim; ++i) {
      const double d = deriv[i] * force_constant * conf_it->second;
      DEBUG(1,"\tderiv[" << i+i << "]: " << d);
      coordinates[i]->apply(d);
    }
  }
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
        break;
    }
    cutoff_rel[i] = conversion_factor_grid * cutoff[i];
    width_rel[i] = conversion_factor_grid * width[i];
    grid_min_rel[i] = conversion_factor * grid_min[i];
    grid_max_rel[i] = conversion_factor * grid_max[i];

    // check for weird periodicity settings
    if (grid_max_rel[i] != grid_min_rel[i])
      grid_spacing_rel[i] = num_grid_points[i] / (grid_max_rel[i] - grid_min_rel[i]);
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
  DEBUG(1, "Calculating coordinates of umbrella " << id);
  for (unsigned int i = 0; i < dim(); ++i) {
    coordinates[i]->calculate(conf);
  }
}

std::string util::Umbrella::str() const {
  std::ostringstream os;
  os << "\tUmbrella (" << id << "): ";
  if (enabled) {
    os << "enabled, " << (building ? "building" : "frozen") << ".";
  } else {
    os << "disabled.";
  }
  os.precision(8);
  os << "\n\tforce constant: " << std::setw(15) << force_constant;
  os << "\n\tNumber of attached coordinates: " << dim() << "\n";
  for(unsigned int i = 0; i < dim(); ++i) {
    os << "\t\t- " << coordinates[i]->str() << "\n";
    os << "\t\t    function    : ";
    switch(functional_form[i]) {
      case ff_gaussian :
        os << "gaussian";
        break;
      case ff_polynomial :
        os << "truncated polynomial";
        break;
      default:
        os << "unkown";
    }
    os << "\n";
    os << "\t\t    width       : " << std::setw(15) << width[i] << "\n"
       << "\t\t    cutoff      : " << std::setw(15) << cutoff[i] << "\n"
       << "\t\t    grid points : " << num_grid_points[i] << "\n"
       << "\t\t    grid min    : " << std::setw(15) << grid_min[i] << "\n"
       << "\t\t    grid max    : " << std::setw(15) << grid_max[i] << "\n"
       << "\t\t    " << configurations.size() << " configurations visited.\n";
  }

  return os.str();
}

