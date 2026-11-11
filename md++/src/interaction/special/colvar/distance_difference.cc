/** Geometry and derivatives of a signed distance difference. */
#include "../../../stdheader.h"
#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../math/periodicity.h"
#include "../../../util/template_split.h"
#include "distance_difference.h"

namespace {

template<math::boundary_enum B>
void calculate_distance_difference(topology::Topology &topo,
    configuration::Configuration &conf, interaction::Distance_Difference_Colvar &cv,
    int &error)
{
  math::Periodicity<B> periodicity(conf.current().box);
  cv.ct = 0.0;
  cv.strain_derivative = 0.0;
  std::fill(cv.derivatives.begin(), cv.derivatives.end(), math::Vec(0));

  for (size_t pair = 0; pair < 2; ++pair) {
    const size_t first = 2 * pair;
    math::Vec displacement;
    periodicity.nearest_image(cv.atoms[first]->pos(conf, topo),
                              cv.atoms[first + 1]->pos(conf, topo), displacement);
    const double distance = math::abs(displacement);
    if (distance <= 1e-12) {
      io::messages.add("DISTANCE_DIFFERENCE: zero-length distance vector.",
                       "Distance_Difference_Colvar", io::message::error);
      return;
    }
    const double sign = pair == 0 ? 1.0 : -1.0;
    cv.ct += sign * distance;
    cv.derivatives[first] = sign * displacement / distance;
    cv.derivatives[first + 1] = -cv.derivatives[first];
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        cv.strain_derivative(a, b) += displacement(a) * cv.derivatives[first](b);
      }
    }
  }
  error = 0;
}

} // anonymous namespace

int interaction::Distance_Difference_Colvar::init(topology::Topology &topo,
    configuration::Configuration &conf, simulation::Simulation &sim,
    std::ostream &os, bool quiet)
{
  if (params == NULL || params->endpoints.size() != 4) {
    io::messages.add("DISTANCE_DIFFERENCE requires four endpoints.",
                     "Distance_Difference_Colvar", io::message::error);
    return 1;
  }
  atoms.clear();
  for (auto &endpoint : params->endpoints) atoms.push_back(&endpoint);
  derivatives.resize(4);
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));
  if (!quiet) os << "Distance-difference colvar restraint interaction\n";
  return 0;
}

int interaction::Distance_Difference_Colvar::calculate(topology::Topology &topo,
    configuration::Configuration &conf, simulation::Simulation &sim)
{
  int error = 1;
  SPLIT_BOUNDARY(calculate_distance_difference, topo, conf, *this, error);
  return error;
}
