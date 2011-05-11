/**
 * @file le_coordinate.cc
 * collective coordinates for LE-US
 */

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../math/periodicity.h"

#include "../util/le_coordinate.h"
#include "../util/template_split.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE leus

util::LE_Coordinate * util::LE_Dihedral_Coordinate::
clone() const {
  return new LE_Dihedral_Coordinate(umbrella_id(), i,j,k,l);
}

void util::LE_Dihedral_Coordinate::
calculate(configuration::Configuration & conf) {
  SPLIT_BOUNDARY(_calculate, conf);
}

template<math::boundary_enum B>
void util::LE_Dihedral_Coordinate::
_calculate(configuration::Configuration & conf) {
  // create the periodicity
  // this is maybe slow but may be improved later if required
  math::Periodicity<B> periodicity(conf.current().box);

  // save the configuration pointer
  m_conf = &conf;

  // calculate the dihedral angle and it's derivative by the coordinates
  math::VArray & pos = conf.current().pos;
  math::Vec rij, rkj, rkl, rmj, rnk, rim, rln;

  periodicity.nearest_image(pos(i), pos(j), rij);
  periodicity.nearest_image(pos(k), pos(j), rkj);
  periodicity.nearest_image(pos(k), pos(l), rkl);

  //calculate phi, cross- and dot-products
  rmj = cross(rij, rkj);
  rnk = cross(rkj, rkl);
  const double dkj2 = abs2(rkj);
  const double dmj2 = abs2(rmj);
  const double dnk2 = abs2(rnk);
  
  assert(dkj2 != 0.0);
  const double frim = dot(rij, rkj) / dkj2;
  const double frln = dot(rkl, rkj) / dkj2;

  rim = rij - frim * rkj;
  rln = frln * rkj - rkl;
  const double dim = sqrt(abs2(rim));
  const double dln = sqrt(abs2(rln));

  const double ip = dot(rim, rln);
  double cosphi = ip / (dim * dln);

  if (cosphi < -1.0) cosphi = -1.0;
  if (cosphi > 1.0) cosphi = 1.0;

  // get the value
  phi = acos(cosphi);

  const double sign = dot(rij, rnk);

  if (sign < 0) {
    phi = -phi;
  }

  // get the derivative
  assert(dmj2 != 0.0);
  assert(dnk2 != 0.0);
  const math::Vec dphi_dri = (sqrt(dkj2) / dmj2) * rmj;
  const math::Vec dphi_drl = -(sqrt(dkj2) / dnk2) * rnk;
  const math::Vec dphi_drj = (frim - 1) * dphi_dri - frln*dphi_drl;
  const math::Vec dphi_drk = -1.0 * dphi_dri - dphi_drj - dphi_drl;

  // save the derivates
  fi = - dphi_dri;
  fj = - dphi_drj;
  fk = - dphi_drk;
  fl = - dphi_drl;
}

double util::LE_Dihedral_Coordinate::
get_value(double grid_min, double grid_max) const {
  double angle = phi;
  DEBUG(10, "angle: " << phi << " grid: " << grid_min << "-" << grid_max);
  // map it on the grid
  while (angle < grid_min) angle += 2.0 * math::Pi;
  while (angle >= grid_max) angle -= 2.0 * math::Pi;
  return angle;
}

double util::LE_Dihedral_Coordinate::
get_deviation(const double & grid_value) const {
  double diff = phi - grid_value;
  // map it on the grid
  while (diff >= math::Pi) diff -= 2.0 * math::Pi;
  while (diff < -math::Pi) diff += 2.0 * math::Pi;
  return diff;
}

void util::LE_Dihedral_Coordinate::
apply(double deriv) {
  // just multiply the derivative with the negative derivative of the
  // angle by the atom positions to get the force.
  const math::Vec & force_i = deriv * fi;
  DEBUG(10, "force(" << i+1 << "): " << math::v2s(force_i));
  m_conf->current().force(i) += force_i;
  const math::Vec & force_j = deriv * fj;
  DEBUG(10, "force(" << j+1 << "): " << math::v2s(force_j));
  m_conf->current().force(j) += force_j;
  const math::Vec & force_k = deriv * fk;
  DEBUG(10, "force(" << k+1 << "): " << math::v2s(force_k));
  m_conf->current().force(k) += force_k;
  const math::Vec & force_l= deriv * fl;
  DEBUG(10, "force(" << l+1 << "): " << math::v2s(force_l));
  m_conf->current().force(l) += force_l;
}

std::string util::LE_Dihedral_Coordinate::
str() const {
  std::ostringstream os;
  os << "LE Dihedral (" << i+1 << "-" << j+1 << "-" << k+1 << "-" << l+1 << ")";
  return os.str();
}

