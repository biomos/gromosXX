/**
 * @file coordnum.cc
 * @brief Implementation of coordination number collective variable
 */

#include <limits>
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../math/periodicity.h"

#include "../../interaction/interaction_types.h"
#include "../../interaction/special/colvar/colvar.h"
#include "../../interaction/special/colvar/coordnum.h"
#include "../../interaction/special/colvar/colvar_utils.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

using namespace interaction;

template<math::boundary_enum B, math::virial_enum V>
static int _calculate_coordnum_colvar(
  topology::Topology &topo,
  configuration::Configuration &conf,
  simulation::Simulation &sim,
  math::VArray &derivatives,
  topology::coordnum_restraint_struct *params,
  double &ct)
{
  math::Periodicity<B> periodicity(conf.current().box);
  ct = 0;
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  for (size_t i = 0; i < (*params).atoms1.size(); i++) {
    for (size_t j = 0; j < (*params).atoms2.size(); j++) {
      math::Vec v;
      double dfunc, func;
      periodicity.nearest_image(
        (*params).atoms1[i].pos(conf, topo),
        (*params).atoms2[j].pos(conf, topo),
        v);
      double rdist = math::abs(v);

      func = switching_function(rdist / (*params).rcut, dfunc, (*params).nn, (*params).mm);
      ct += func;

      math::Vec d = dfunc * v / rdist;
      derivatives[i] += d;
      derivatives[(*params).atoms1.size() + j] -= d;
    }
  }
  return 0;
}

int Coordnum_Colvar::calculate(topology::Topology &topo,
                               configuration::Configuration &conf,
                               simulation::Simulation &sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_coordnum_colvar,
    topo, conf, sim, derivatives, params, ct);
  return 0;
}

int Coordnum_Colvar::init(topology::Topology &topo,
                          configuration::Configuration &conf,
                          simulation::Simulation &sim,
                          std::ostream &os,
                          bool quiet)
{
  targetvalue = (*params).cont0;
  rcut = (*params).rcut;
  mm = (*params).mm;
  nn = (*params).nn;
  w0 = (*params).w0;

  for (auto &a : (*params).atoms1)
    atoms.push_back(&a);
  for (auto &a : (*params).atoms2)
    atoms.push_back(&a);

  derivatives.resize((*params).atoms1.size() + (*params).atoms2.size());

  if (!quiet)
    os << "Coordnum restraint interaction\n";
  return 0;
}
