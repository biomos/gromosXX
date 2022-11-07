#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../interaction/interaction.h"

#include "../util/debug.h"

#include "bs_potentials.h"
#include "bs_vector.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE bs_leus

void util::BS_Potential::calcPotAndForce(double gridPointCoordinate,
        double &potential,
        double &force) {
  m_activeGridPoint = int (gridPointCoordinate + 0.5); // NINT(Gamma / R * r)
  int gridPointIndex1 = int(gridPointCoordinate);
  int gridPointIndex2 = int(gridPointCoordinate) + 1;
  double x = gridPointCoordinate - gridPointIndex1;
  double x2 = x * x;
  double x3 = x * x2;
  potential = m_memory[gridPointIndex1] * (1 - 3 * x2 + 2 * x3);
  potential += m_memory[gridPointIndex2] * (3 * x2 - 2 * x3);
  force = (m_memory[gridPointIndex1] - m_memory[gridPointIndex2]) * 6 * (x2 - x);
  DEBUG(10, "Grid Point Coordinate: " << gridPointCoordinate << "; active Grid Point " << m_activeGridPoint);
  DEBUG(10, "index1: " << gridPointIndex1 << "; index2: " << gridPointIndex2);
  DEBUG(10, "Potential: " << potential << "; force: " << force);
}

double util::BS_Potential::calcWeight(std::vector<double> &potentials, double beta) {
  double sum = 0;
  std::vector<double>::iterator it = potentials.begin(),
          to = potentials.end();
  for (; it != to; it++) {
    sum += exp(-beta * (*it - m_potential));
  }
  m_weight = 1.0 / sum;
  return m_weight;
}

bool util::BS_Potential::updateMemory(bool updateAuxMem) {
  DEBUG(12, "I will " << updateAuxMem << " update auxmem");
  DEBUG(10, "I update the memory by " << m_weight << " * " << m_forceIncrement);
  m_memory[m_activeGridPoint] += m_weight * m_forceIncrement;

  if (updateAuxMem) {
    DEBUG(10, "Update auxiliary memory by " << m_weight);
    m_auxiliaryMemory[m_activeGridPoint] += m_weight;

    // are all memory points bigger than local visiting cutoff?
    bool bigger = true;
    std::vector<double>::iterator it = m_auxiliaryMemory.begin(),
            to = m_auxiliaryMemory.end();
    for (; it != to; it++) {
      bigger = bigger && (*it >= m_localCutoff);
      if (!bigger)
        return false;
    }
    return bigger; // should be true
  }
  return false;
}

void util::BS_Potential::setMemoryToZero() {
  m_memory.assign(num_gp, 0);
}

void util::BS_Potential::setAuxMemoryToZero() {
  m_auxiliaryMemory.assign(num_gp, 0);
}

void util::BS_Potential::setMemory(std::vector<double>& newMemory) {
  if (newMemory.size() != num_gp) {
    std::ostringstream os;
    os << "Cannot set memory to new memory. Size of new Memory ("
            << newMemory.size() << ") is not equal the number of grid points ("
            << num_gp << ") in potential with id " << id << "!";
    io::messages.add(os.str(), "BS_Potential", io::message::error);
    return;
  }
  m_memory = newMemory;
}

void util::BS_Potential::setAuxMemory(std::vector<double>& newMemory) {
  if (newMemory.size() != num_gp) {
    std::ostringstream os;
    os << "Cannot set memory to new memory. Size of new Auxiliary Memory ("
            << newMemory.size() << ") is not equal the number of grid points ("
            << num_gp << ") in potential with id " << id << "!";
    io::messages.add(os.str(), "BS_Potential", io::message::error);
    return;
  }
  m_auxiliaryMemory = newMemory;
}

void util::BS_Potential::getAuxMemory(std::vector<double>& newMemory) {
  newMemory = m_auxiliaryMemory;
}

void
util::BS_Potential::setMemoryParameters(double forceIncrement, int localCutoff) {
  m_forceIncrement = forceIncrement;
  m_localCutoff = localCutoff;
}

std::string util::BS_Potential::traj_str() {
  std::ostringstream os;
  os << std::setw(3) << id << " "
          << std::setw(3) << m_potentialType
          << std::setw(12) << m_potential << " "
          << std::setw(12) << m_weight << " "
          << std::setw(4) << m_activeGridPoint;

  return os.str();
}

/********************************************************
  THE SPHERE
 */
double util::BS_Sphere::calcPotential(BS_Vector & bs_pos) {
  // vec(r_k)
  DEBUG(8, "\nSphere " << id << " Center: " << m_center.str() << " Position: " << bs_pos.str());
  BS_Vector radialVector;
  bs_pos.minus(m_center, radialVector);
  double radialDistance = radialVector.normalize();
  DEBUG(6, "; RadialVector |" << radialVector.str() << "| = " << radialDistance);

  if (radialDistance >= m_radius) {
    double diff2 = radialDistance - m_radius;
    m_potDerivatives = radialVector * (diff2 * m_force_const);
    DEBUG(10, "Multiply radial Vector with diff^2 * force = " << (diff2 * m_force_const));
    diff2 = diff2 * diff2;
    m_potential = m_memory[num_gp - 1];
    m_potential += m_half_force_const * diff2;
    m_activeGridPoint = (num_gp - 1);
    DEBUG(8, "Outside: Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Sphere " << id << " Derivatives = " << m_potDerivatives.str());
  } else {
    double gridPointCoordinate = m_gridPointScaling * radialDistance;
    //double gridPointCoordinate = (num_gp - 1) * pow((radialDistance / m_radius), radialVector.size());
    DEBUG(6, "GridPointCoordinate = " << gridPointCoordinate);
    double force = 0.0;
    calcPotAndForce(gridPointCoordinate, m_potential, force);
    DEBUG(10, "Force = " << force);
    m_potDerivatives = radialVector * (m_gridPointScaling * force);
    //m_potDerivatives = radialVector * (gridPointCoordinate * radialVector.size() / radialDistance * force);
    DEBUG(6, "Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Sphere " << id << " Derivatives = " << m_potDerivatives.str());
  }
  return m_potential;
}

std::string util::BS_Sphere::str() {
  std::ostringstream os;
  os << std::setw(3) << id
          << std::setw(10) << m_force_const << " "
          << std::setw(6) << m_radius << " "
          << m_center.str();
  return os.str();
}

/*************************************************
 THE STICK
 */
double util::BS_Stick::calcPotential(BS_Vector& bs_pos) {
  BS_Vector relativePosition;
  bs_pos.minus(m_startPoint, relativePosition);
  double longitudinalDistance = relativePosition.dot(m_unitLongitudinal);
  DEBUG(8, "\nStick " << id << "; relativePosition " << relativePosition.str() << "longDist " << longitudinalDistance);
  DEBUG(8, "\tUnit Longitudinal Vector: " << m_unitLongitudinal.str());
  // At the beginning of the stick
  if (longitudinalDistance <= 0) {
    m_distance = relativePosition.normalize();
    m_activeGridPoint = 0;
    m_potential = m_memory[0];
    double offset = m_distance - m_half_width;
    DEBUG(8, "offset = " << offset);
    if (offset > 0) {
      m_potential += m_half_force_const * offset * offset;
      m_potDerivatives = relativePosition * (offset * m_force_const);
    }
    else {
      m_potDerivatives.assign(relativePosition.size(), 0);
    }
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Stick " << id << " Derivatives = " << m_potDerivatives.str());
    //return false;
  }    // At the end of the stick
  else if (longitudinalDistance >= m_length) {
    bs_pos.minus(m_endPoint, relativePosition);
    m_distance = relativePosition.normalize();
    m_activeGridPoint = num_gp - 1;
    m_potential = m_memory[num_gp - 1];
    double offset = m_distance - m_half_width;
    DEBUG(8, "offset = " << offset);
    if (offset > 0) {
      m_potential += m_half_force_const * offset * offset;
      m_potDerivatives = relativePosition * (offset * m_force_const);
    }
    else {
      m_potDerivatives.assign(relativePosition.size(), 0);
    }
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Stick " << id << " Derivatives = " << m_potDerivatives.str());
  }
    // aside of the stick
  else {
    double gridPointCoordinate = m_gridPointScaling * longitudinalDistance;
    double force = 0.0;
    calcPotAndForce(gridPointCoordinate, m_potential, force);
    m_potDerivatives = m_unitLongitudinal * (force * m_gridPointScaling);
    BS_Vector transversalVector;
    relativePosition.minus(m_unitLongitudinal * longitudinalDistance,
            transversalVector);
    m_distance = transversalVector.normalize();
    DEBUG(8, "Transversal distance: " << m_distance << "; from : " << transversalVector.str());
    double offset = m_distance - m_half_width;
    if (offset > 0) {
      m_potential += m_half_force_const * offset * offset;
      m_potDerivatives += transversalVector * (offset * m_force_const);
    }
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Stick " << id << " Derivatives = " << m_potDerivatives.str());
  }
  return m_potential;
}

std::string util::BS_Stick::str() {
  std::ostringstream os;
  os << std::setw(3) << id
          << std::setw(10) << m_force_const << " "
          << std::setw(6) << m_half_width << " "
          << m_startPoint.str() << " "
          << m_endPoint.str();
  return os.str();
}

/*************************************************
 THE SNAKE
 */
util::BS_Snake::BS_Snake(int id, int num_gp, double force_const,
        std::vector<util::BS_Vector> points, double halfwidth) :
util::BS_Potential(id, (points.size() - 1) * (num_gp + 1) + 1, force_const) {
  int num_sticks = points.size() - 1;
  DEBUG(5, "Snake " << id << " with " << num_sticks << " sticks.");
  for (int i = 0; i < num_sticks; i++) {
    int where = util::BS_Distorted_Stick::in_between;
    if (i == 0) {
      DEBUG(6, "Adding Distorted stick at start (" << where << ")");
      where = util::BS_Distorted_Stick::start;
    } else if (i == (num_sticks - 1)) {
      where = util::BS_Distorted_Stick::end;
      DEBUG(6, "Adding Distorted stick at end (" << where << ")");
    } else {
      DEBUG(6, "Adding Distorted stick in between (" << where << ")");
    }
    m_sticks.push_back(util::BS_Distorted_Stick(i, num_gp + 2, force_const,
            points[i], points[i + 1], halfwidth, where));
  }
  for (int i = 0; i < num_sticks; i++) {
    int prev_i = i -1;
    int next_i = i + 1;
    if (i == 0){
      prev_i++;
    } else if (i == (num_sticks - 1)){
      next_i--;
    }
    m_sticks[i].init_vectors(m_sticks[prev_i].direction(), m_sticks[next_i].direction());
  }
  m_potentialType = bs_snake;
}

double util::BS_Snake::calcPotential(BS_Vector& bs_pos) {
  DEBUG(5, "Calculate the Potential of snake " << id);
  int size = m_sticks.size();
  int mem_stride = (num_gp - 1) / size;
  std::vector<double>::iterator mem = m_memory.begin();
  double distances[size];
  bool interacts[size];

  double min_distance = 0.0;
  int min_i = -1;

  for (int i = 0; i < size; i++) {
    std::vector<double> stick_mem(mem, mem + mem_stride + 1);
    m_sticks[i].setMemory(stick_mem);
    mem = mem_stride + mem;
    interacts[i] = m_sticks[i].interacts(bs_pos);
    if (interacts[i]){
      distances[i] = m_sticks[i].distance();
      min_distance = distances[i];
      DEBUG(4, distances[i] << " from stick " << i);
    } else {
      DEBUG(4, "No interaction with stick " << i);

    }
  }

  bool any_interacting = false;

  for (int i = 0; i < size; i++){
    DEBUG(7, "Does it interact with stick " << i << "? " << (interacts[i] ? "yes" : "no"));
    any_interacting = any_interacting || interacts[i];
    if (interacts[i]){
      min_distance = distances[i];
      min_i = i;
    }
  }

  for (int i = 0; i < size; i++) {
    if (interacts[i] && distances[i] < min_distance) {
      min_distance = distances[i];
      min_i = i;
    }
  }
  DEBUG(4, "Minimal distance is from stick " << min_i);

  if (!any_interacting || min_i == -1){
    DEBUG(4, "No interacting stick found!");
    std::ostringstream os;
    os << "Couldn't find a stick interacting with snake ";
    os << id << std::endl;
    io::messages.add(os.str(), "BS_Potential", io::message::error);
    return 0;
  }
    
  m_potential = m_sticks[min_i].calcPotential(bs_pos);
  m_sticks[min_i].getDerivatives(m_potDerivatives);
  m_activeGridPoint = min_i * mem_stride + m_sticks[min_i].getActiveGridPoint();
  // find internal gp
  return m_potential;
}

std::string util::BS_Snake::str() {
  std::ostringstream os;
  os << "# SNAKE  ID  SIZE\n    " << std::setw(3) 
          << id << "   " << m_sticks.size() << "\n";
  for (unsigned int i = 0; i < m_sticks.size(); i++) {
    os << "\t" << m_sticks[i].str() << "\n";
  }
  os << "# END SNAKE";
  return os.str();
}

/*************************************************
 THE DISTORTED STICK
 */

util::BS_Distorted_Stick::BS_Distorted_Stick(int id, int num_gp, double force_const,
        BS_Vector startPoint, BS_Vector endPoint,
        double half_width, int is_where) :
    util::BS_Stick(id, num_gp, force_const, startPoint, endPoint, half_width) 
{
  switch (is_where) {
    case 0: {
      where = in_between;
      DEBUG(6, "Create stick in between. (" << id << ")");
      break;
    }
    case 1: {
      where = start;
       DEBUG(6, "Create stick at start. (" << id << ")");
     break;
    }
    case 2: {
      where = end;
      DEBUG(6, "Create stick at end. (" << id << ")");
      break;
    }
    default: {
      std::ostringstream os;
      os << "Wrong value for 'where' in Distorted Stick (" << where << "). ";
      os << "Only 0, 1 or 2." << std::endl;
      io::messages.add(os.str(), "BS_Potential", io::message::error);
    }
  }
}
    
    
void
util::BS_Distorted_Stick::init_vectors(BS_Vector prev_dir, BS_Vector next_dir) {
  DEBUG(5, "Initilaize vectors for distorted stick " << id);
  BS_Vector S, T, Sigma, Tau;
  double start_tangens = 0.0, end_tangens = 0.0;
  // K = prev_dir
  // L = m_unitLongitudinal
  // M = next_dir

  if (where != start) {
    // S = (K - L) / |K - L|
    S = prev_dir - m_unitLongitudinal;
    double lengthS = S.normalize();
    // T = (K - (K * L)L) / |K - (K * L)L|
    T = prev_dir - m_unitLongitudinal * m_unitLongitudinal.dot(prev_dir);
    T.normalize();
    if (lengthS == 0.0){
      S = T;
    }
    // (S*L) / (S*T)
    start_tangens = S.dot(m_unitLongitudinal);
    if (start_tangens != 0.0){
      start_tangens /= S.dot(T);
    }
    // N_S = L - (L * S)L  # start_norm
    //m_unitLongitudinal.minus(S * m_unitLongitudinal.dot(S), m_start_norm);
    //m_start_norm.normalize();
  }

  if (where == end) {
    Tau = T;
    Sigma = Tau;
    end_tangens = 0.0;
  } else {
    // Sigma (L - M) / |L -M|
    Sigma = m_unitLongitudinal - next_dir;
    double lengthSigma = Sigma.normalize();
    // Tau = (M - (M*L)L) / |M - (M*L)L|
    Tau = next_dir - m_unitLongitudinal * m_unitLongitudinal.dot(next_dir);
    Tau.normalize();
    if (lengthSigma == 0.0){
      Sigma = Tau;
    }
    //  (Sigma * L) / (Sigma * Tau)
    end_tangens = Sigma.dot(m_unitLongitudinal);
    if (end_tangens != 0.0){
      end_tangens /= Sigma.dot(Tau);
    }
    // N_E = - (L - (L * Sigma)L)  # end_norm
    //m_unitLongitudinal.minus(S * m_unitLongitudinal.dot(Sigma), m_end_norm);
    //m_end_norm.scale(-1);
    //m_end_norm.normalize();
  }
  
  if (where == start) {
    T = Tau;
    S = T;
    start_tangens = 0.0;
  } 

  DEBUG(7, "L = " << m_unitLongitudinal.str());
  DEBUG(7, "K = " << prev_dir.str());
  DEBUG(7, "M = " << next_dir.str());
  DEBUG(7, "S = " << S.str());
  DEBUG(7, "T = " << T.str());
  DEBUG(7, "Sigma = " << Sigma.str());
  DEBUG(7, "Tau = " << Tau.str());
  DEBUG(7, "SL / ST = " << start_tangens);
  DEBUG(7, "Sigma L / Sigma Tau = " << end_tangens);  
  
  // h = r * h_base
  // H = U + r * H_base
  // h_base = L - T * (S*L) / (S*T)
  m_h_base = m_unitLongitudinal - T * start_tangens;
  // H_base = Tau * (Sigma * L) / (Sigma * Tau) - T * (S*L) / (S*T)
  m_H_base = Tau * end_tangens - T * start_tangens;
  
  DEBUG(7, "h_base = " << m_h_base.str());
  DEBUG(7, "H_base = " << m_H_base.str());
 
}

bool util::BS_Distorted_Stick::interacts(BS_Vector& bs_pos) {
  DEBUG(5, "Does system interact with distorted stick " << id << "?");
  BS_Vector r; // relative position from the stick start
  bs_pos.minus(m_startPoint, r);
  BS_Vector rho; // relative position from the stick start
  bs_pos.minus(m_endPoint, rho);
  DEBUG(6, "Sizes of vectors: ");
  DEBUG(6, "r: " << r.size());
  DEBUG(6, "rho: " << rho.size());
  DEBUG(6, "unit L: " << m_unitLongitudinal.size());
  /**if ((where == start && m_end_norm.dot(rho) < 0) ||
          (where == end && m_start_norm.dot(r) < 0) ||
          (m_end_norm.dot(rho) < 0 && m_start_norm.dot(r) < 0)) {
    return false;
  }*/
  m_at_end = true;
  double longitudinalDistance = m_unitLongitudinal.dot(r);
  DEBUG(7, "Longitudinal Distance: " << longitudinalDistance);
  if (where == start && longitudinalDistance < 0) {
    m_distance = r.normalize();
    radialVec = r;
    //return true;
  } else if (where == end && longitudinalDistance >= m_length) {
    m_distance = rho.normalize();
    radialVec = rho;
    //return true;
  } else {
    m_at_end = false;
    BS_Vector perpDistance;
    r.minus(m_unitLongitudinal * longitudinalDistance, perpDistance);

    m_h = r.dot(m_h_base);
    m_H = m_length + r.dot(m_H_base);
    m_h_over_H = m_h / m_H;
    
    DEBUG(8, "r = " << r.str());
    DEBUG(7, "h = " << m_h << "; H = " << m_H);
    DEBUG(7, "h / H: " << m_h_over_H);
    if (m_h_over_H < 0.0 || m_h_over_H > 1.0){
      return false;
    }
    
    m_deriv_h_over_H = m_h_base - (m_H_base * m_h_over_H);
    m_deriv_h_over_H.scale(1.0 / m_H);

    //r.minus(m_unitLongitudinal * (m_h_over_H * m_length), m_unit_tilted_perp);
    //m_tilted_perp_dist = m_unit_tilted_perp.normalize();
    //m_distance = m_tilted_perp_dist;
    r.minus(m_unitLongitudinal * longitudinalDistance, m_unit_perp);
    m_perp_dist = m_unit_perp.normalize();
    m_distance = m_perp_dist;
  }
  return true;
}

double util::BS_Distorted_Stick::calcPotential(BS_Vector& bs_pos) {
  DEBUG(5, "Calculate the potential of distorted stick " << id);
  if (m_at_end) {
    double offset = m_distance - m_half_width;
    double memory = 0.0;
    if (where == start) {
      memory = m_memory[0];
    } else if (where == end) {
      memory = m_memory[num_gp - 1];
    } else {
      std::ostringstream os;
      os << "It looks like I am not at the ends of a distorted stick, "
              << "but still have to calculate the potential for it. "
              << "ID is " << id << std::endl;
      io::messages.add(os.str(), "BS_Potential", io::message::error);
    }
    if (offset > 0) {
      m_potential = m_half_force_const * offset * offset + memory;
      m_potDerivatives = radialVec * (offset * m_force_const);
    } else {
      m_potential = memory;
      m_potDerivatives.assign(bs_pos.size(), 0);
    }
  } else { // in between
    double gridPointCoordinate = (num_gp - 1)* m_h_over_H;
    double force = 0.0;
    calcPotAndForce(gridPointCoordinate, m_potential, force);

    m_potDerivatives = m_deriv_h_over_H * (force * (num_gp - 1));
    DEBUG(5, "Distorted Stick Derivative from grid: " << m_potDerivatives.str());

    double offset = m_distance - m_half_width;
    if (offset > 0) {
      m_potential += m_half_force_const * offset * offset;

      //m_tilted_perp_deriv = m_unit_tilted_perp * (1 - m_length / m_H);
      //m_potDerivatives += m_tilted_perp_deriv * (m_force_const * offset);
      m_potDerivatives += m_unit_perp * (m_force_const * offset);
    }
  }

  DEBUG(5, "Total Distorted Stick Derivative: " << m_potDerivatives.str());
  return m_potential;
}

std::string util::BS_Distorted_Stick::str() {
  std::ostringstream os;
  os << "Distorted Stick ";
  if (where == start){
    os << "at start.   ";
  } else if (where == end){
    os << "at end.     ";
  } else {
    os << "in between. ";
  }
  os << std::setw(3) << id
          << std::setw(10) << m_force_const << " "
          << std::setw(6) << m_half_width << " "
          << m_startPoint.str() << " "
          << m_endPoint.str();
  return os.str();
}

/*************************************************
 THE DISTORTED STICK
 */

util::BS_Pipe::BS_Pipe(int id, int num_gp_long, int num_gp_perp, double force_const,
            util::BS_Pipe_Param start, util::BS_Pipe_Param end) :
      BS_Potential(id, num_gp_long * num_gp_perp, force_const),
      m_start(start), m_end(end),
      m_num_long_gp(num_gp_long), m_num_perp_gp(num_gp_perp)
{
  m_unitLongitudinal = m_end.point - m_start.point;
  m_length = m_unitLongitudinal.normalize();
  
  m_inner_slope = (m_end.inner_width - m_start.inner_width) / m_length;
  m_outer_slope = (m_end.outer_width - m_start.outer_width) / m_length;
  DEBUG(5, "Slopes: inner: " << m_inner_slope << "; outer: " << m_outer_slope);
  DEBUG(5, "grid points " << num_gp_long << " x " << num_gp_perp);
  
  m_long_conversion = (num_gp_long - 1) / m_length;
  m_potentialType = bs_pipe;
}

double util::BS_Pipe::calcPotential(BS_Vector &bs_pos){
  BS_Vector rel_pos = bs_pos - m_start.point;
  DEBUG(8, "Position: " << bs_pos.str());
  DEBUG(8, "Relative Position: " << rel_pos.str());
  
  // LONGITUDINAL
  double long_dist = rel_pos.dot(m_unitLongitudinal);
  BS_Vector long_pos = m_unitLongitudinal * long_dist;
  
  double long_pot = 0.0;
  m_potDerivatives.assign(bs_pos.size(), 0);
  where_in_potential longitudinal;  
  double l_diff = long_dist;
  if (l_diff < 0){
    longitudinal = below_lower;
    long_pot = l_diff * l_diff * m_half_force_const;
    m_potDerivatives += m_unitLongitudinal * (m_force_const * l_diff);
  } else if ((l_diff = long_dist - m_length) > 0){
    longitudinal = above_upper;
    long_pot = l_diff * l_diff * m_half_force_const;
    m_potDerivatives += m_unitLongitudinal * (m_force_const * l_diff);
  } else {
    longitudinal = inside;
  }  
  DEBUG(8, "long dist = " << long_dist);
  DEBUG(8, "Longitudinal: " << longitudinal);
  DEBUG(8, "long pot: " << long_pot);
  
  // PERPENDICULAR
  BS_Vector perp_pos = rel_pos - long_pos;
  double perp_dist = perp_pos.normalize();
  
  double inner_cutoff = m_inner_slope * long_dist + m_start.inner_width;
  double outer_cutoff = m_outer_slope * long_dist + m_start.outer_width;
  
  double perp_pot = 0.0;
  where_in_potential perpendicular;
  double inner_p_diff = perp_dist - inner_cutoff;
  double outer_p_diff = perp_dist - outer_cutoff;
  if (inner_p_diff < 0){
    perpendicular = below_lower;
    perp_pot = inner_p_diff * inner_p_diff * m_half_force_const;
    m_potDerivatives += (perp_pos - m_unitLongitudinal * m_inner_slope) * 
            (m_force_const * inner_p_diff);
  } else if (outer_p_diff > 0){
    perpendicular = above_upper;
    perp_pot = outer_p_diff * outer_p_diff * m_half_force_const;
    m_potDerivatives += (perp_pos - m_unitLongitudinal * m_outer_slope) * 
                            (m_force_const * outer_p_diff);
  } else {
    perpendicular = inside;
  }
  DEBUG(8, "perp dist = " << perp_dist);
  DEBUG(8, "Perpendicular: " << perpendicular);
  DEBUG(8, "perp pot: " << perp_pot);

  // GRID
  double long_grid = m_long_conversion * long_dist;
  int i = 0;
  int active_i = 0;
  if (longitudinal == above_upper){
    i = m_num_long_gp - 1;
    active_i = i;
  } else if (longitudinal == inside) {
    i = int(long_grid);
    active_i = int(long_grid + 0.5);
  }
  
  const double cutoff_diff = outer_cutoff - inner_cutoff;
  const double perp_conversion = (m_num_perp_gp - 1) / cutoff_diff;
  double perp_grid = perp_conversion * (perp_dist - inner_cutoff);
  int j = 0;
  int active_j = 0;
  if (perpendicular == above_upper){
    j = m_num_perp_gp - 1;
    active_j = j;
  } else if (perpendicular == inside){
    j = int(perp_grid);
    active_j = int(perp_grid + 0.5);
  }
  m_activeGridPoint = active_j * m_num_long_gp + active_i;
  DEBUG(8, "Grid Coordinate (l,p): (" << long_grid << ", " << perp_grid << ")");
  DEBUG(8, "active: i: " << active_i << " j: " << active_j << " overall: " << m_activeGridPoint);

  double deltas[4] = {0, 0, 0, 0};
  double ddeltas_di[4] = {0, 0, 0, 0}; // derivative of delta_i times delta_j
  double ddeltas_dj[4] = {0, 0, 0, 0}; // derivative of delta_j times delta_i
  if (longitudinal == inside){
    double d_i = long_grid - i;
    double d_i2 = d_i * d_i;
    double delta_1 = 3.0 * d_i2 + 2.0 * d_i2 * d_i;
    double delta_0 = 1.0 - delta_1;
    
    deltas[0] =  delta_0;
    deltas[1] = delta_0;
    deltas[2] = delta_1;
    deltas[3] = delta_1;
    
    double ddelta = 6 * (d_i - d_i2);
    ddeltas_di[0] = -ddelta;
    ddeltas_di[1] = -ddelta;
    ddeltas_di[2] =  ddelta;
    ddeltas_di[3] =  ddelta;
  }
  else if (longitudinal == above_upper) {
    deltas[2] = 1.0;
    deltas[3] = 1.0;
  } else { // below
    deltas[0] = 1.0;
    deltas[1] = 1.0;
  }
  
  for (int i = 0; i < 4; i++){
    ddeltas_dj[i] = deltas[i];
  }   
  
  if (perpendicular == inside){
    double d_j = perp_grid - j;
    double d_j2 = d_j * d_j;
    double delta_1 = 3.0 * d_j2 + 2.0 * d_j2 * d_j;
    double delta_0 = 1.0 - delta_1;   
    
    deltas[0] *= delta_0;
    deltas[2] *= delta_0;
    deltas[1] *= delta_1;
    deltas[3] *= delta_1; 
    
    double ddelta = 6 * (d_j - d_j2);
    ddeltas_dj[0] *= -ddelta;
    ddeltas_dj[2] *= -ddelta;
    ddeltas_dj[1] *=  ddelta;
    ddeltas_dj[3] *=  ddelta;
            
    ddeltas_di[0] *= delta_0;
    ddeltas_di[2] *= delta_0;
    ddeltas_di[1] *= delta_1;
    ddeltas_di[3] *= delta_1; 
  }
  else if (perpendicular == above_upper) {
    deltas[0] = 0.0;
    deltas[2] = 0.0;
    ddeltas_di[0] = 0.0;
    ddeltas_di[2] = 0.0;
  } else {
    deltas[1] = 0.0;
    deltas[3] = 0.0;
    ddeltas_di[1] = 0.0;
    ddeltas_di[3] = 0.0;
  }
  
  bool has_delta[4];
  int index = j * m_num_long_gp + i;
  double grid_pot = 0.0;
  for (int m = 0; m < 2; m++){    // m <- j
    for (int n = 0; n < 2; n++){  // n <- i
      DEBUG(8, "detla_" << n << "," << m << ": " << deltas[m * 2 + n])
      if ((has_delta[m * 2 + n] = (deltas[m * 2 + n]) != 0)){
        grid_pot += deltas[m * 2 + n] * m_memory[index + m * m_num_long_gp + n];
      }
    }
  }
  DEBUG (8, "Grid Potential = " << grid_pot);

  BS_Vector grid_deriv;
  grid_deriv.assign(bs_pos.size(), 0.0);
  if (longitudinal == inside) {
    double pot_deriv = 0.0;
    for (int m = 0; m < 2; m++) {   // m <- j
      for (int n = 0; n < 2; n++) { // n <- i
        DEBUG(8, "d detla_" << n << "," << m << "/ di: " << ddeltas_di[m * 2 + n])
        if (has_delta[m * 2 + n]) {
          pot_deriv += ddeltas_di[m * 2 + n] * m_memory[index + m * m_num_long_gp + n];
        }
      }
    }
    if (pot_deriv != 0.0){
      grid_deriv += m_unitLongitudinal * (pot_deriv * m_long_conversion);
    }
  }

  if (perpendicular == inside) {
    double pot_deriv = 0.0;
    for (int m = 0; m < 2; m++) {   // m <- j
      for (int n = 0; n < 2; n++) { // n <- i
        DEBUG(8, "d detla_" << n << "," << m << "/ dj: " << ddeltas_dj[m * 2 + n])
        if (has_delta[m * 2 + n]) {
          pot_deriv += ddeltas_dj[m * 2 + n] * m_memory[index + m * m_num_long_gp + n];
        }
      }
    }
    if (pot_deriv != 0.0){
      double d_j_prefactor = perp_conversion / cutoff_diff;
      d_j_prefactor *= (m_inner_slope * outer_p_diff - m_outer_slope * inner_p_diff);
      grid_deriv += m_unitLongitudinal * d_j_prefactor + perp_pos * perp_conversion;
    }
  }
  
  m_potential = grid_pot + perp_pot + long_pot;
  m_potDerivatives += grid_deriv;
  DEBUG (8, "Total Potential = " << m_potential);
  
  return m_potential;
}

std::string util::BS_Pipe::str() {
  std::ostringstream os;
  os << std::setw(3) << id
          << std::setw(10) << m_force_const << " "
          << std::setw(6) << m_start.inner_width << " "
          << std::setw(6) << m_start.outer_width << " "
          << std::setw(6) << m_end.inner_width << " "
          << std::setw(6) << m_end.outer_width << " "
          << std::setw(6) << m_start.point.str() << " "
          << std::setw(6) << m_end.point.str() << " ";
  return os.str();
}
