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


#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE bs_leus

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
  const double lowerLimit = 0.1;
  const double upperLimit = 6.1;
  DEBUG(5, "Close to Zero = " << closeToZero << "; Close To Periode = " << closeToPeriode);
  DEBUG(5, "phi = " << phi);
  if (closeToZero && (phi > math::Pi)){
    phi -= 2 * math::Pi;
    DEBUG(5, "Decreasing phi...");
  }
  else if (closeToPeriode && (phi < math::Pi)){
    DEBUG(5, "Increasing phi...");
    phi += 2 * math::Pi;
  }
  DEBUG(5, "after: phi = " << phi);
  if (phi < lowerLimit){
    closeToZero = true;
    closeToPeriode = false;
    DEBUG(5, "I am below the lower limit");
  } else if (phi > upperLimit) {
    closeToZero = false;
    closeToPeriode = true;
    DEBUG(5, "I am over the upper Limit");
  } else {
    DEBUG(5, "I am in between");
    closeToZero = false;
    closeToPeriode = false;
  }
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
  double deriv = derivatives.begin()->value * m_rad2degree;
  DEBUG(10, "addForce: deriv: " << deriv << "; orig. value: " << derivatives.begin()->value);
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

void util::BS_Dihedral::getInternalCoordinates(std::vector<BS_Dimension>& coord) const {
  coord.clear();
  BS_Dimension dihedralCoord;
  dihedralCoord.value = red_phi * m_rad2degree;
  dihedralCoord.periodicity = m_periodicity;
  coord.push_back(dihedralCoord);
}

std::string util::BS_Dihedral::str() const {
  std::ostringstream os;
  os << "Dihedral: (" << i + 1 << ", " << j + 1 << ", " << k + 1 << ", "
          << j + 1 << ")\n";
  return os.str();
}
/**************************************************
 * BS_VECTOR
 */
double util::BS_Vector::abs2(){
  double length2 = 0;
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++){
    length2 += it->value * it->value;
  }
  return length2;
}

double util::BS_Vector::normalize(){
  double length = sqrt(this->abs2());
  double inverseLength = 1 / length;
  this->scale(inverseLength);
  return length;
}

void util::BS_Vector::nullify(){
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++){
    it->value = 0;
  }
}

void util::BS_Vector::scale(const double scalar){
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++){
    DEBUG(10, "Value * scalar =  " << it->value << " * " << scalar);
    it->value = it->value * scalar;
    DEBUG(10, "\t = " << it->value)
  }
}

void util::BS_Vector::minus(const BS_Vector& subtrahend, BS_Vector& result){
  if (this->size() != subtrahend.size()){
    io::messages.add("Two BS_Vectors with different sizes (minus)!", "BS_Vector",
            io::message::critical);
    return;
  }
  result.clear();
  BS_Vector::iterator m_i = this->begin(),
          m_end = this->end();
  BS_Vector::const_iterator s_i = subtrahend.begin();
  
  BS_Dimension diff;
  for (; m_i != m_end; m_i++, s_i++){
    diff.value = m_i->value - s_i->value;
    DEBUG(12, "BS_Vector.minus: Value: " << diff.value);
    diff.periodicity = 0.0;
    if(m_i->periodicity != 0.0){
      DEBUG(12, "BS_Vector.minus: Periodicity: " << m_i->periodicity);
      diff.periodicity = m_i->periodicity;
      double halfPeriod = 0.5 * m_i->periodicity;
      while (diff.value > halfPeriod)
        diff.value -= m_i->periodicity;
      while (diff.value < -halfPeriod)
        diff.value += m_i->periodicity;
    }
    DEBUG(12, "BS_Vector.minus: After applying periodicity: " << diff.value);
    result.push_back(diff);
  }
}

util::BS_Vector util::BS_Vector::operator *(const double scalar){
  BS_Vector result;
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  BS_Dimension scaled;
  for (; it != to; it++){
    scaled.value = it->value * scalar;
    scaled.periodicity = it->periodicity;
    // Not checking for periodicity: ok?
    result.push_back(scaled);
  }
  return result;
}

util::BS_Vector util::BS_Vector::operator +(const BS_Vector& summand) {
  if (this->size() != summand.size()){
    io::messages.add("Two BS_Vectors with different sizes (plus)!", "BS_Vector",
            io::message::critical);
  }
  BS_Dimension sum;
  BS_Vector result;
  BS_Vector::const_iterator s_i = summand.begin();
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++, s_i++){
    sum.value  = it->value + s_i->value;
    sum.periodicity = it->periodicity;
    if (sum.periodicity != 0.0){
      while (sum.value > sum.periodicity)
        sum.value -= sum.periodicity;
      while (sum.value < 0)
        sum.value += sum.periodicity;
    }
    result.push_back(sum);
  }
  return result;
}

void util::BS_Vector::operator +=(const BS_Vector& summand) {
  if (this->size() != summand.size()){
    io::messages.add("Two BS_Vectors with different sizes (+=)!", "BS_Vector",
            io::message::critical);
  }
  BS_Vector::const_iterator s_i = summand.begin();
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++, s_i++){
    it->value += s_i->value;
    if (it->periodicity != 0.0){
      while (it->value > it->periodicity)
        it->value -= it->periodicity;
      while (it->value < 0)
        it->value += it->periodicity;
    }
  }
}

double util::BS_Vector::dot(const BS_Vector &other){
  if (this->size() != other.size()){
    io::messages.add("Two BS_Vectors with different sizes (dot)!", "BS_Vector",
            io::message::critical);
  }
  double dotProduct = 0;
  BS_Vector::const_iterator other_i = other.begin();
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++, other_i++){
    dotProduct += it->value * other_i->value;
  }
  return dotProduct;
}

void util::BS_Vector::create(std::vector<double>& values, 
        std::vector<double>& periodicities)
{
  this->clear();
  util::BS_Dimension dim;
  std::vector<double>::iterator val_i = values.begin(),
          val_end = values.end(),
          period_i = periodicities.begin();
  for (; val_i != val_end; val_i++, period_i++){
    dim.value = *val_i;
    dim.periodicity = *period_i;
    this->push_back(dim);
  }
}

std::string util::BS_Vector::str(){
  std::ostringstream os;
  os << "(";
  BS_Vector::const_iterator vec_i = this->begin(),
          vec_end = this->end();
  
  if (vec_i != vec_end){
    os << vec_i->value;
    vec_i++;
  }
  else {
    os << "'empty'";
  }
  for (; vec_i != vec_end; vec_i++){
    os << ", " << vec_i->value;
  }
  os << ")";
  return os.str();
}
