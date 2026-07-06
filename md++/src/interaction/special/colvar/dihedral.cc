/**
 * @file dihedral.cc
 * @brief Implementation of dihedral collective variable
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
#include "../../interaction/special/colvar/dihedral.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

using namespace interaction;

template<math::boundary_enum B, math::virial_enum V>
static int _calculate_dihedral_colvar(
  topology::Topology &topo,
  configuration::Configuration &conf,
  simulation::Simulation &sim,
  math::VArray &derivatives,
  topology::dihedral_restraint_struct *params,
  double &ct)
{
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  math::Vec rij, rkj, rkl;
  periodicity.nearest_image(pos(params->i), pos(params->j), rij);
  periodicity.nearest_image(pos(params->k), pos(params->j), rkj);
  periodicity.nearest_image(pos(params->k), pos(params->l), rkl);

  const math::Vec rmj = cross(rij, rkj);
  const math::Vec rnk = cross(rkj, rkl);

  const double dkj2 = abs2(rkj);
  const double dmj2 = abs2(rmj);
  const double dnk2 = abs2(rnk);
  const double dkj = std::sqrt(dkj2);
  const double dmj = std::sqrt(dmj2);
  const double dnk = std::sqrt(dnk2);

  if (dkj <= 1.0e-12 || dmj <= 1.0e-12 || dnk <= 1.0e-12) {
    io::messages.add("DIHEDRAL colvar: derivative is singular.",
                     "Dihedral_Colvar", io::message::warning);
    ct = 0.0;
    return 0;
  }

  double acs = dot(rmj, rnk) / (dmj * dnk);
  if (acs > 1.0) acs = 1.0;
  if (acs < -1.0) acs = -1.0;

  ct = std::acos(acs);
  if (dot(rij, rnk) < 0.0) {
    ct *= -1.0;
  }

  const math::Vec di = dkj / dmj2 * rmj;
  const math::Vec dl = -dkj / dnk2 * rnk;
  const double dj1 = dot(rij, rkj) / dkj2 - 1.0;
  const double dj2 = dot(rkl, rkj) / dkj2;

  derivatives[0] = di;
  derivatives[3] = dl;
  derivatives[1] = dj1 * derivatives[0] - dj2 * derivatives[3];
  derivatives[2] = -1.0 * (derivatives[0] + derivatives[1] + derivatives[3]);

  return 0;
}

int Dihedral_Colvar::calculate(topology::Topology &topo,
                               configuration::Configuration &conf,
                               simulation::Simulation &sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_dihedral_colvar,
    topo, conf, sim, derivatives, params, ct);
  return 0;
}

int Dihedral_Colvar::init(topology::Topology &topo,
                          configuration::Configuration &conf,
                          simulation::Simulation &sim,
                          std::ostream &os,
                          bool quiet)
{
  targetvalue = params->phi;
  w0 = params->w0;

  m_atoms.clear();
  m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->i)));
  m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->j)));
  m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->k)));
  m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->l)));

  atoms.clear();
  for (size_t i = 0; i < m_atoms.size(); ++i) {
    atoms.push_back(&m_atoms[i]);
  }

  derivatives.resize(4);
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  if (!quiet)
    os << "Dihedral colvar restraint interaction\n";
  return 0;
}
