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
//-------------------------------------
// LE_Distance_coordinate
//-------------------------------------

util::LE_Coordinate * util::LE_Distance_Coordinate::
clone() const {
  return new LE_Distance_Coordinate(umbrella_id(), i,j);
}

void util::LE_Distance_Coordinate::
calculate(configuration::Configuration & conf) {
  SPLIT_BOUNDARY(_calculate, conf);
}

template<math::boundary_enum B>
void util::LE_Distance_Coordinate::
_calculate(configuration::Configuration & conf) {
  // create the periodicity
  // this is maybe slow but may be improved later if required
  math::Periodicity<B> periodicity(conf.current().box);

  // save the configuration pointer
  m_conf = &conf;
  // calculate the distance and it's derivative by the coordinates
  math::VArray & pos = conf.current().pos;
  math::Vec rij;

  periodicity.nearest_image(pos(i), pos(j), rij);

  // get the value
  dist = abs(rij);
//  dist = rij.abs();


  // save the derivates
  fi = math::Vec(1.0);
  fj = math::Vec(-1.0);
}

double util::LE_Distance_Coordinate::
get_value(double grid_min, double grid_max) const {
  double d = dist;
  DEBUG(10, "distance: " << dist << " grid: " << grid_min << "-" << grid_max);
  // map it on the grid
  // while (angle < grid_min) angle += 2.0 * math::Pi;
  // while (angle >= grid_max) angle -= 2.0 * math::Pi;
  // what to do for a distance that is beyond the grid?
  return d;
}

double util::LE_Distance_Coordinate::
get_deviation(const double & grid_value) const {
  double diff = dist - grid_value;
  // map it on the grid
  // check if it is outside the grid??
  return diff;
}

void util::LE_Distance_Coordinate::
apply(double deriv) {
  // just multiply the derivative with the negative derivative of the
  // distance by the atom positions to get the force.
  const math::Vec & force_i = deriv * fi;
  DEBUG(10, "force(" << i+1 << "): " << math::v2s(force_i));
  m_conf->current().force(i) += force_i;
  const math::Vec & force_j = deriv * fj;
  DEBUG(10, "force(" << j+1 << "): " << math::v2s(force_j));
  m_conf->current().force(j) += force_j;
}

std::string util::LE_Distance_Coordinate::
str() const {
  std::ostringstream os;
  os << "LE Distance (" << i+1 << "-" << j+1 << ")";
  return os.str();
}

//-------------------------------------
// LE_DistanceField_coordinate
//-------------------------------------

// we need some helper functions
int neighbour(int i, int j, std::vector<int> &ngrid);
double assignment_1d(const int &p, const double &xi);

util::LE_DistanceField_Coordinate::LE_DistanceField_Coordinate(int id, topology::Topology &topo, simulation::Simulation &sim){

  // take the (virtual) atoms from the distancefield.
  va_i = topo.disfield_restraints().v1;
  va_j = topo.disfield_restraints().v2;

  m_topo = &topo;
  m_sim = &sim;

  umbrella_id(id);
  
}

util::LE_Coordinate * util::LE_DistanceField_Coordinate::
clone() const {
  return new LE_DistanceField_Coordinate(umbrella_id(), *m_topo, *m_sim);
}

void util::LE_DistanceField_Coordinate::
calculate(configuration::Configuration & conf) {
  SPLIT_BOUNDARY(_calculate, conf);
}

template<math::boundary_enum B>
void util::LE_DistanceField_Coordinate::
_calculate(configuration::Configuration & conf) {

  // store the configuration for later
  m_conf = &conf;
  
  // this is mostly a copy of the _calculate_distance_field_interaction function
  //math::Vec f;
  std::vector<int> &ngrid = conf.special().distancefield.ngrid;
  math::Box &box = conf.current().box;
  double grid = m_sim->param().distancefield.grid;
  std::vector<double> &distance = conf.special().distancefield.distance;
  
  math::Vec pos_j = va_j.pos(conf, *m_topo);
  //pos(sim.param().distancefield.atom_j);
  
  // determine the nearest grid center
  // first we get the position of the particle in terms of grid coordinates (double)
  // we also get from that the coordinates of the lowest gridpoint of the eight surrounding points;
  math::Vec gpos_j;
  math::GenericVec<int> grid_j;
  
  for(int i=0; i<3; i++){
    gpos_j[i] = (pos_j[i] + box(i,i)/2)/grid;
    grid_j[i] = int(gpos_j[i]);
  }

  // fill an array with the indices of the eight neighbouring gridpoints
  std::vector<int> eightpoints(8);
  eightpoints[0] = grid_j[2] * ngrid[0] * ngrid[1] + grid_j[1] * ngrid[0] + grid_j[0];
  eightpoints[1] = neighbour(eightpoints[0], 1, ngrid);
  eightpoints[2] = neighbour(eightpoints[0], 3, ngrid);
  eightpoints[3] = neighbour(eightpoints[1], 3, ngrid);
  eightpoints[4] = neighbour(eightpoints[0], 5, ngrid);
  eightpoints[5] = neighbour(eightpoints[1], 5, ngrid);
  eightpoints[6] = neighbour(eightpoints[2], 5, ngrid);
  eightpoints[7] = neighbour(eightpoints[3], 5, ngrid);

  DEBUG(9, "DF CALC, eightpoints\n\t" << eightpoints[0] << " " << distance[eightpoints[0]] << "\n\t"
	<< eightpoints[1] << " " << distance[eightpoints[1]] << "\n\t"
	<< eightpoints[2] << " " << distance[eightpoints[2]] << "\n\t"
	<< eightpoints[3] << " " << distance[eightpoints[3]] << "\n\t"
	<< eightpoints[4] << " " << distance[eightpoints[4]] << "\n\t"
	<< eightpoints[5] << " " << distance[eightpoints[5]] << "\n\t"
	<< eightpoints[6] << " " << distance[eightpoints[6]] << "\n\t"
	<< eightpoints[7] << " " << distance[eightpoints[7]]);

  // and one with their coordinates in grids
  std::vector<math::Vec> eightpos(8);
  eightpos[0][0] = double(grid_j[0]);
  eightpos[0][1] = double(grid_j[1]);
  eightpos[0][2] = double(grid_j[2]);
  eightpos[1] = eightpos[0] + math::Vec(1.0,0.0,0.0);
  eightpos[2] = eightpos[0] + math::Vec(0.0,1.0,0.0);
  eightpos[3] = eightpos[0] + math::Vec(1.0,1.0,0.0);
  eightpos[4] = eightpos[0] + math::Vec(0.0,0.0,1.0);
  eightpos[5] = eightpos[0] + math::Vec(1.0,0.0,1.0);
  eightpos[6] = eightpos[0] + math::Vec(0.0,1.0,1.0);
  eightpos[7] = eightpos[0] + math::Vec(1.0,1.0,1.0);

  // calculate the derivatives in the point. For this, we need to go even beyond
  // the points. Hmmm. Does this carry the risk of getting a force from the protein
  // at very long distances? But it has to be if we want to get a smooth surface.
  std::vector<math::Vec> eightderiv(8);
  for(unsigned int i=0; i< eightderiv.size(); i++){
    eightderiv[i][0] = (distance[neighbour(eightpoints[i], 1, ngrid)] - 
			distance[neighbour(eightpoints[i], 0, ngrid)]);
    eightderiv[i][1] = (distance[neighbour(eightpoints[i], 3, ngrid)] - 
			distance[neighbour(eightpoints[i], 2, ngrid)]); 
    eightderiv[i][2] = (distance[neighbour(eightpoints[i], 5, ngrid)] - 
			distance[neighbour(eightpoints[i], 4, ngrid)]); 
    eightderiv[i] /= (2*grid);
    DEBUG(9, "DF CALC, 8 derivative: " << i << " : " << eightderiv[i][0] << " " << eightderiv[i][1] << " " << eightderiv[i][2]);
  }
  
  // assign the average distance using an assignment function of order 2.
  dist=0;
  math::Vec deriv(0.0,0.0,0.0);
  
  for(unsigned int i=0; i< eightpoints.size(); i++){
    double P = assignment_1d(2,gpos_j[0] - eightpos[i][0]) * 
      assignment_1d(2,gpos_j[1] - eightpos[i][1]) *
      assignment_1d(2,gpos_j[2] - eightpos[i][2]);
    dist += distance[eightpoints[i]] * P;
    deriv += eightderiv[i] * P;
  }

  // and this is where we are done. We calculated the distance and the derivative of the 
  // field with respect to the coordinates, which we store here in the vectors fi and fj

  fi= deriv;
  fj=-deriv;
}

double util::LE_DistanceField_Coordinate::
get_value(double grid_min, double grid_max) const {
  double d = dist;
  DEBUG(10, "distance: " << dist << " grid: " << grid_min << "-" << grid_max);
  // map it on the grid
  // while (angle < grid_min) angle += 2.0 * math::Pi;
  // while (angle >= grid_max) angle -= 2.0 * math::Pi;
  // what to do for a distance that is beyond the grid?
  return d;
}

double util::LE_DistanceField_Coordinate::
get_deviation(const double & grid_value) const {
  double diff = dist - grid_value;
  // map it on the grid
  // check if it is outside the grid??
  return diff;
}

void util::LE_DistanceField_Coordinate::
apply(double deriv) {
  // just multiply the derivative with the negative derivative of the
  // distance by the atom positions to get the force.

  const math::Vec & force_i = deriv * fi;
  DEBUG(10, "force(i): " << math::v2s(force_i));
  va_i.force(*m_conf, *m_topo, force_i);
  //m_conf->current().force(i) += force_i;

  const math::Vec & force_j = deriv * fj;
  DEBUG(10, "force(j): " << math::v2s(force_j));
  //m_conf->current().force(j) += force_j;
  va_j.force(*m_conf, *m_topo, force_j);
  
}

std::string util::LE_DistanceField_Coordinate::
str() const {
  std::ostringstream os;
  os << "LE DistanceField according to DISTANCEFIELD input";
  return os.str();
}
