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
  const std::vector<util::Virtual_Atom*> &atoms,
  double &ct,
  int &error,
  math::Matrix &strain_derivative)
{
  error = 0;
  strain_derivative = 0.0;
  math::Periodicity<B> periodicity(conf.current().box);
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  math::Vec rij, rkj, rkl;
  periodicity.nearest_image(atoms[0]->pos(conf,topo), atoms[1]->pos(conf,topo), rij);
  periodicity.nearest_image(atoms[2]->pos(conf,topo), atoms[1]->pos(conf,topo), rkj);
  periodicity.nearest_image(atoms[2]->pos(conf,topo), atoms[3]->pos(conf,topo), rkl);

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

  // Reconstruct l relative to j along the same j-k-l image chain.
  const math::Vec rlj = rkj - rkl;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      strain_derivative(a, b) = rij(a) * derivatives[0](b)
                              + rkj(a) * derivatives[2](b)
                              + rlj(a) * derivatives[3](b);
    }
  }

  return 0;
}

int Dihedral_Colvar::calculate(topology::Topology &topo,
                               configuration::Configuration &conf,
                               simulation::Simulation &sim)
{
  int error = 1; // Also fail if boundary/virial dispatch rejects the setup.
  SPLIT_VIRIAL_BOUNDARY(_calculate_dihedral_colvar,
    topo, conf, sim, derivatives, atoms, ct, error, strain_derivative);
  return error;
}

int Dihedral_Colvar::init(topology::Topology &topo,
                          configuration::Configuration &conf,
                          simulation::Simulation &sim,
                          std::ostream &os,
                          bool quiet)
{
  targetvalue = params->phi;
  w0 = params->w0;

  if (endpoints && endpoints->size()!=4) return 1;
  m_atoms.clear();
  if (endpoints) m_atoms=*endpoints;
  else {
    m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->i)));
    m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->j)));
    m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->k)));
    m_atoms.push_back(util::Virtual_Atom(util::va_explicit, std::vector<int>(1, params->l)));
  }

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
