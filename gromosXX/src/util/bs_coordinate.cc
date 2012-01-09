#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../interaction/interaction.h"

#include "../math/periodicity.h"
#include "../io/message.h"

#include "../util/template_split.h"
#include "../util/debug.h"

#include "bs_coordinate.h"
#include "bs_vector.h"


#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE bs_leus

//=============== DIHEDRAL ===================================
void util::BS_Dihedral::
calculateInternalCoord(configuration::Configuration & conf) {
  DEBUG(8, "Calculate Internal Coordinates according to pbc");
  SPLIT_BOUNDARY(_calculate, conf);
}

template<math::boundary_enum B>
void util::BS_Dihedral::
_calculate(configuration::Configuration & conf) {
  // This is more or less the same code as in LEUS
  DEBUG(8, "Dihedral: Calculate the Internal Coordinates of Dihedral angle");
  // create the periodicity
  // this is maybe slow but may be improved later if required
  math::Periodicity<B> periodicity(conf.current().box);

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
    // After phi == pi, phi decreases again, thus:
    phi = -phi + 2 * math::Pi; // Put it into the range of [0; 360]
  }
  double diff = old_phi - phi;
  if (diff > math::Pi){ // we just moved from just below 360 to just above 0
    m_counter++;
  } 
  else if (diff < - math::Pi){ // vice versa
    m_counter--;
  }
  old_phi = phi;
  phi += m_counter * 2 * math::Pi;
  red_phi = phi / m_red_fac;  

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

void util::BS_Dihedral::addForces(configuration::Configuration& conf, 
        BS_Vector& derivatives){
  assert(derivatives.size() == m_dimension);
  DEBUG(10, "Derivative: " << derivatives.str());
  double deriv = derivatives[0] * m_rad2degree / m_red_fac;
  DEBUG(10, "addForce: deriv: " << deriv << "; orig. value: " << derivatives[0]);
  // just multiply the derivative with the negative derivative of the
  // angle by the atom positions to get the force.
  const math::Vec & force_i = deriv * fi;
  DEBUG(10, "force(" << i+1 << "): " << math::v2s(force_i));
  conf.current().force(i) += force_i;
  DEBUG(10, "force_total(" << i+1 << "): " << math::v2s(conf.current().force(i)));
  const math::Vec & force_j = deriv * fj;
  DEBUG(10, "force(" << j+1 << "): " << math::v2s(force_j));
  conf.current().force(j) += force_j;
  const math::Vec & force_k = deriv * fk;
  DEBUG(10, "force(" << k+1 << "): " << math::v2s(force_k));
  conf.current().force(k) += force_k;
  const math::Vec & force_l= deriv * fl;
  DEBUG(10, "force(" << l+1 << "): " << math::v2s(force_l));
  conf.current().force(l) += force_l;
}

void util::BS_Dihedral::getInternalCoordinates(BS_Vector& coord) const {
  coord.clear();
  coord.push_back(red_phi * m_rad2degree);
}

std::string util::BS_Dihedral::str() const {
  std::ostringstream os;
  os << "Dihedral: (" << i + 1 << ", " << j + 1 << ", " << k + 1 << ", "
          << j + 1 << ")\n";
  return os.str();
}

// =================== DISTANCE ====================
void util::BS_Distance::
calculateInternalCoord(configuration::Configuration & conf) {
  DEBUG(8, "Calculate Internal Coordinates according to pbc");
  SPLIT_BOUNDARY(_calculate, conf);
}

template<math::boundary_enum B>
void util::BS_Distance::
_calculate(configuration::Configuration & conf) {
  DEBUG(8, "Dihedral: Calculate the Internal Coordinates of Distance");
    // create the periodicity
  // this is maybe slow but may be improved later if required
  math::Periodicity<B> periodicity(conf.current().box);

  // calculate the dihedral angle and it's derivative by the coordinates
  math::VArray & pos = conf.current().pos;
  math::Vec rij;

  periodicity.nearest_image(pos(i), pos(j), rij);
  
  m_distance = abs(rij);
  m_reducedDistance = m_distance / m_red_fac;
  
  fi = rij / m_distance;
  fj = - fi;
}

void util::BS_Distance::addForces(configuration::Configuration& conf, 
        BS_Vector& derivatives){
  assert(derivatives.size() == m_dimension);
  DEBUG(10, "Derivative: " << derivatives.str());
  //double deriv = derivatives.begin()->value;
  double deriv = derivatives[0];
  DEBUG(10, "addForce: derivative: " << deriv);
  // just multiply the derivative with the negative derivative of the
  // angle by the atom positions to get the force.
  const math::Vec & force_i = deriv * fi;
  DEBUG(10, "force(" << i+1 << "): " << math::v2s(force_i));
  conf.current().force(i) += force_i;
  DEBUG(10, "force_total(" << i+1 << "): " << math::v2s(conf.current().force(i)));
  const math::Vec & force_j = deriv * fj;
  DEBUG(10, "force(" << j+1 << "): " << math::v2s(force_j));
  conf.current().force(j) += force_j;
}

void util::BS_Distance::getInternalCoordinates(BS_Vector& coord) const {
  coord.clear();
  coord.push_back(m_reducedDistance);
}

std::string util::BS_Distance::str() const {
  std::ostringstream os;
  os << "Distance: (" << i + 1 << ", " << j + 1 << ")\n";
  return os.str();
}