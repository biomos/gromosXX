/**
 * @file distance.cc
 * @brief Implementation of distance collective variable
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../math/periodicity.h"

#include "../../interaction/interaction_types.h"
#include "../../interaction/special/colvar/colvar.h"
#include "../../interaction/special/colvar/distance.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

using namespace interaction;

static int _distance_rah_mask(int rah, math::Vec &mask, int &rah_mode)
{
  if (rah >= -1 && rah <= 1) {
    mask = math::Vec(1, 1, 1);
    rah_mode = rah;
    return 0;
  }
  if (rah >= 9 && rah <= 11) {
    mask = math::Vec(1, 1, 0);
    rah_mode = rah - 10;
    return 0;
  }
  if (rah >= 19 && rah <= 21) {
    mask = math::Vec(1, 0, 1);
    rah_mode = rah - 20;
    return 0;
  }
  if (rah >= 29 && rah <= 31) {
    mask = math::Vec(0, 1, 1);
    rah_mode = rah - 30;
    return 0;
  }
  if (rah >= 39 && rah <= 41) {
    mask = math::Vec(1, 0, 0);
    rah_mode = rah - 40;
    return 0;
  }
  if (rah >= 49 && rah <= 51) {
    mask = math::Vec(0, 1, 0);
    rah_mode = rah - 50;
    return 0;
  }
  if (rah >= 59 && rah <= 61) {
    mask = math::Vec(0, 0, 1);
    rah_mode = rah - 60;
    return 0;
  }

  return 1;
}

template<math::boundary_enum B, math::virial_enum V>
static int _calculate_distance_colvar(
  topology::Topology &topo,
  configuration::Configuration &conf,
  simulation::Simulation &sim,
  math::VArray &derivatives,
  topology::distance_restraint_struct *params,
  const math::Vec &dim_mask,
  double &ct)
{
  math::Periodicity<B> periodicity(conf.current().box);

  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  math::Vec v;
  periodicity.nearest_image(params->v1.pos(conf, topo),
                            params->v2.pos(conf, topo),
                            v);

  for (int i = 0; i < 3; ++i)
    v[i] *= dim_mask[i];

  ct = math::abs(v);

  if (ct > 1.0e-12) {
    derivatives[0] =  v / ct;
    derivatives[1] = -v / ct;
  }

  return 0;
}

int Distance_Colvar::calculate(topology::Topology &topo,
                               configuration::Configuration &conf,
                               simulation::Simulation &sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_distance_colvar,
    topo, conf, sim, derivatives, params, dim_mask, ct);
  return 0;
}

int Distance_Colvar::init(topology::Topology &topo,
                          configuration::Configuration &conf,
                          simulation::Simulation &sim,
                          std::ostream &os,
                          bool quiet)
{
  targetvalue = params->r0;
  w0 = params->w0;

  atoms.clear();
  atoms.push_back(&params->v1);
  atoms.push_back(&params->v2);

  derivatives.resize(2);
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  if (_distance_rah_mask(params->rah, dim_mask, rah_mode)) {
    std::ostringstream msg;
    msg << "Do not know how to handle " << params->rah
        << " for RAH in distance colvar" << std::endl;
    io::messages.add(msg.str(), "Distance_Colvar", io::message::error);
    return 1;
  }

  if (!quiet)
    os << "Distance colvar restraint interaction\n";

  return 0;
}