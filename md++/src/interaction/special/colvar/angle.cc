/**
 * @file angle.cc
 * @brief Implementation of angle collective variable
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
#include "../../interaction/special/colvar/angle.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

using namespace interaction;

template<math::boundary_enum B, math::virial_enum V>
static int _calculate_angle_colvar(
  topology::Topology &topo,
  configuration::Configuration &conf,
  simulation::Simulation &sim,
  math::VArray &derivatives,
  topology::angle_restraint_struct *params,
  double &ct)
{
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  math::Vec rij, rkj;
  periodicity.nearest_image(pos(params->i), pos(params->j), rij);
  periodicity.nearest_image(pos(params->k), pos(params->j), rkj);

  const double dij = math::abs(rij);
  const double dkj = math::abs(rkj);
  if (dij <= 1.0e-12 || dkj <= 1.0e-12) {
    io::messages.add("ANGLE colvar: zero-length angle vector.",
                     "Angle_Colvar", io::message::error);
    ct = 0.0;
    return 1;
  }

  double cost = dot(rij, rkj) / (dij * dkj);
  if (cost > 1.0) cost = 1.0;
  if (cost < -1.0) cost = -1.0;

  ct = std::acos(cost);

  const double sint = std::sqrt(std::max(0.0, 1.0 - cost * cost));
  if (sint <= 1.0e-12) {
    io::messages.add("ANGLE colvar: derivative is singular near 0 or pi.",
                     "Angle_Colvar", io::message::warning);
    return 0;
  }

  const math::Vec dcost_di = (rkj / dkj - rij / dij * cost) / dij;
  const math::Vec dcost_dk = (rij / dij - rkj / dkj * cost) / dkj;

  derivatives[0] = -dcost_di / sint;
  derivatives[2] = -dcost_dk / sint;
  derivatives[1] = -derivatives[0] - derivatives[2];

  return 0;
}

int Angle_Colvar::calculate(topology::Topology &topo,
                            configuration::Configuration &conf,
                            simulation::Simulation &sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_angle_colvar,
    topo, conf, sim, derivatives, params, ct);
  return 0;
}

int Angle_Colvar::init(topology::Topology &topo,
                       configuration::Configuration &conf,
                       simulation::Simulation &sim,
                       std::ostream &os,
                       bool quiet)
{
  targetvalue = params->theta;
  w0 = params->w0;

  m_atoms.clear();
  m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->i)));
  m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->j)));
  m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->k)));

  atoms.clear();
  for (size_t i = 0; i < m_atoms.size(); ++i) {
    atoms.push_back(&m_atoms[i]);
  }

  derivatives.resize(3);
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  if (!quiet)
    os << "Angle colvar restraint interaction\n";
  return 0;
}
