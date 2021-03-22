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

// ============== BASE: COORDINATE =========================
std::string
util::BS_Coordinate::init_str()
{
  std::ostringstream os;
  os << std::setw(3) << m_id << " "
          << std::setw(8) << m_red_fac << " "
          << std::setw(5) << m_dimension << " "
          << std::setw(7) << m_type;
  return os.str();
}

//=============== DIHEDRAL ===================================
util::BS_Dihedral::BS_Dihedral(int id, unsigned int i, unsigned int j,
            unsigned int k, unsigned int l, double red_fac) :
    i(i), j(j), k(k), l(l) 
{ 
  cid(id);
  m_type = dihedral;
  m_dimension = 1;
  m_red_fac = red_fac;
  m_red_pi = math::Pi / m_red_fac;
  m_red_2pi = 2 * m_red_pi;
  m_rad2degree = 180 / math::Pi;
  m_counter = 0;
  old_phi = math::Pi;
  DEBUG(10, "i: " << i << " j: " << j << " k: " << k << " l: " << l);
}

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
  DEBUG(8, "phi = " << phi << " reduced_phi = " << red_phi);

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

void util::BS_Dihedral::setOldPos(BS_Vector::iterator it, BS_Vector::iterator to)
{
  old_phi = *it * m_red_fac * math::Pi / 180;
  phi = *it * m_red_fac * math::Pi / 180;
  double two_pi = 2 * math::Pi;
  int mod = 0;
  if (old_phi > two_pi){
    do {
      old_phi -= two_pi;
      mod++;
    }
    while (old_phi > two_pi);
  }
  if (old_phi < 0){
    do {
      old_phi += two_pi;
      mod--;
    }
    while (old_phi < 0);
  }
  m_counter = mod;
  DEBUG(8, "SetOldPos: Phi = " << phi << "; counter = " << m_counter);

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
          << l + 1 << ")\n";
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

  // calculate the distances and it's derivative by the coordinates
  math::VArray & pos = conf.current().pos;
  math::Vec rij;

  periodicity.nearest_image(pos(i), pos(j), rij);
  
  m_distance = abs(rij);
  m_reducedDistance = m_distance / m_red_fac;
  
  fj = rij / m_distance;
  fi = - fj;
}

void util::BS_Distance::addForces(configuration::Configuration& conf, 
        BS_Vector& derivatives){
  assert(derivatives.size() == m_dimension);
  DEBUG(10, "Derivative: " << derivatives.str());
  double deriv = derivatives[0];
  DEBUG(10, "addForce: derivative: " << deriv);
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

// ============= DIHEDRAL_SUM ===========================
util::BS_DihedralSum::BS_DihedralSum(int id, 
            unsigned int i, unsigned int j, unsigned int k, unsigned int l, 
            unsigned int ii, unsigned int jj, unsigned int kk, unsigned int ll, 
            double red_fac) :
            i(i), j(j), k(k), l(l), ii(ii), jj(jj), kk(kk), ll(ll), 
            m_phi(0, i, j, k, l, red_fac), m_psi(0, ii, jj, kk, ll, red_fac) 
{
  cid(id);
  m_type = dihedralSum;
  m_red_fac = red_fac;
  m_dimension = 1;
  m_first = true;
  DEBUG(10, "Create coordinate for dihedral sum");
}
                                               
void util::BS_DihedralSum::
calculateInternalCoord(configuration::Configuration& conf)
{
  m_phi.calculateInternalCoord(conf);
  m_psi.calculateInternalCoord(conf);
  m_sum = m_phi.getPhi() + m_psi.getPhi();
  if (m_first){
    m_first = false;
  } else {
    DEBUG(8, "phi + psi = " << m_sum << "; old sum = " << m_old_sum);
    while ((m_sum - m_old_sum) < (-180 / m_red_fac)) {
      m_sum += 360 / m_red_fac;
      m_phi.incrCounter();
    }
    while ((m_sum - m_old_sum) > (180 / m_red_fac)) {
      m_sum -= 360 / m_red_fac;
      m_phi.decrCounter();
    }
  }
  m_old_sum = m_sum;
  DEBUG(8, "Now: phi + psi = " << m_sum);
}

void util::BS_DihedralSum::
getInternalCoordinates(BS_Vector& coord) const 
{
  coord.clear();
  coord.push_back(m_sum);
}

void util::BS_DihedralSum::
setOldPos(BS_Vector::iterator it, BS_Vector::iterator to)
{
  m_old_sum = *it;
  m_first = false;
}

void util::BS_DihedralSum::
addForces(configuration::Configuration& conf, BS_Vector& derivatives)
{
  DEBUG(8, "Add forces to phi");
  m_phi.addForces(conf, derivatives);
  DEBUG(8, "Add forces to psi");
  m_psi.addForces(conf, derivatives);
}

std::string util::BS_DihedralSum::str() const
{
  std::ostringstream os;
  os << "Dihedral Sum: phi = (" << i + 1 << ", " << j + 1 << ", " 
          << k + 1 << ", " << l + 1 << "); psi = (" << ii + 1 << ", " 
          << jj + 1 << ", " << kk + 1 << ", " << ll + 1 << ")\n";
  return os.str();
}

// ============ Cartesian ===================================
void util::BS_Cartesian::
calculateInternalCoord(configuration::Configuration & conf) {
  DEBUG(8, "Calculate Internal Coordinates according to pbc");
  SPLIT_BOUNDARY(_calculate, conf);
  
  /*m_coordinates.clear();
  if (m_allAtoms){
    for (unsigned int i = 0; i < conf.current().pos.size(); i++){
      for (unsigned int j = 0; j < 3; j++)
        m_coordinates.push_back(conf.current().pos(i)[j]);
    }
  }
  else {
    for (unsigned int i = 0; i < m_atoms.size(); i++) {
      for (unsigned int j = 0; j < 3; j++)
        m_coordinates.push_back(conf.current().pos(m_atoms[i])[j]);
    }
  }*/
}

template<math::boundary_enum B>
void util::BS_Cartesian::
_calculate(configuration::Configuration & conf) {
  DEBUG(8, "Calculate the Internal Coordinates of Cartesian");
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray & pos = conf.current().pos;
  math::Vec vec;

  
  m_coordinates.clear();
  if (m_allAtoms){
    for (unsigned int i = 0; i < pos.size(); i++){
      vec = pos(i);
      periodicity.put_into_positive_box(vec);
      for (unsigned int j = 0; j < 3; j++)
        m_coordinates.push_back(vec[j]);
    }
  }
  else {
    for (unsigned int i = 0; i < m_atoms.size(); i++) {
      vec = pos(m_atoms[i]);
      periodicity.put_into_positive_box(vec);
      for (unsigned int j = 0; j < 3; j++)
        m_coordinates.push_back(vec[j]);
    }
  }
}

void util::BS_Cartesian::getInternalCoordinates(BS_Vector& coord) const {
  coord.clear();
  std::vector<double>::const_iterator it = m_coordinates.begin(),
          to = m_coordinates.end();
  coord.insert(coord.begin(), it, to);
  DEBUG(8, "BS_Reference: Added " << coord.size() << " coordinates. Should be " << m_dimension);
  assert(coord.size() == m_dimension);
}

void util::BS_Cartesian::addForces(configuration::Configuration& conf, 
        BS_Vector& derivatives){
  assert(derivatives.size() == m_dimension);
  DEBUG(10, "Reference: Derivative: " << derivatives.str());
  
  BS_Vector::iterator it = derivatives.begin(),
                      to = derivatives.end();
  std::vector<unsigned int>::iterator atom_i = m_atoms.begin();
  for (; it != to; atom_i++){
    math::Vec force_i(0);
    for (int i = 0; i < 3; i++, it++) {
      force_i(i) = -*it;
    }
    DEBUG(10, "force(" << *atom_i + 1 << "): " << math::v2s(force_i));
    conf.current().force(*atom_i) += force_i;
    DEBUG(10, "total force(" << *atom_i + 1 << "): " << math::v2s(conf.current().force(*atom_i)));
    
  }
}

std::string util::BS_Cartesian::str() const {
  std::ostringstream os;
  os << "Cartesian\n";
  return os.str();
}

// ================ BS_Lambda ==============================
void util::BS_Lambda::
calculateInternalCoord(configuration::Configuration & conf) {
  DEBUG(8, "Calculate lambda.");
  
  // TODO: Get the lambda
  // m_lambda  = ... / m_red_fac;
}

void util::BS_Lambda::addForces(configuration::Configuration& conf, 
        BS_Vector& derivatives){
  assert(derivatives.size() == m_dimension);
  DEBUG(10, "Derivative: " << derivatives.str());
  double deriv = derivatives[0];
  
  // TODO: Add biasing force
  // conf.current().lambda_force += - deriv / m_red_fac;
}

void util::BS_Lambda::getInternalCoordinates(BS_Vector& coord) const {
  coord.clear();
  coord.push_back(m_lambda);
}

std::string util::BS_Lambda::str() const {
  std::ostringstream os;
  os << "Lambda\n";
  return os.str();
}
