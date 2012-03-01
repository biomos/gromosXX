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
        double &force){
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

double util::BS_Potential::calcWeight(std::vector<double> &potentials, double beta)
{
  double sum = 0;
  std::vector<double>::iterator it = potentials.begin(),
          to = potentials.end();
  for (; it != to; it++){
    sum += exp(-beta * (*it - m_potential));
  }
  m_weight = 1.0 / sum;
  return m_weight;
}

bool util::BS_Potential::updateMemory(bool updateAuxMem){
  DEBUG(12, "I will " << updateAuxMem << " update auxmem");
  DEBUG(10, "I update the memory by " << m_weight << " * " << m_forceIncrement);
  m_memory[m_activeGridPoint] += m_weight * m_forceIncrement;
  
  if (updateAuxMem){
    DEBUG(10, "Update auxiliary memory by " << m_weight);
    m_auxiliaryMemory[m_activeGridPoint] += m_weight;
    
    // are all memory points bigger than local visiting cutoff?
    bool bigger = true; 
    std::vector<double>::iterator it = m_auxiliaryMemory.begin(),
            to = m_auxiliaryMemory.end();
    for (; it != to; it++){
      bigger = bigger && (*it >= m_localCutoff);
      if (!bigger)
        return false;
    }
    return bigger; // should be true
  }
  return false;
}

void util::BS_Potential::setMemoryToZero(){
  m_memory.assign(num_gp, 0);
}

void util::BS_Potential::setAuxMemoryToZero(){
  m_auxiliaryMemory.assign(num_gp, 0);
}

void util::BS_Potential::setMemory(std::vector<double>& newMemory){
  if (newMemory.size() != num_gp){
    std::ostringstream os;
    os << "Cannot set memory to new memory. Size of new Memory (" 
            << newMemory.size() << ") is not equal the number of grid points ("
            << num_gp << ") in potential with id " << id << "!";
    io::messages.add(os.str(), "BS_Potential", io::message::error);
    return;
  }
  m_memory = newMemory;
}

void util::BS_Potential::setAuxMemory(std::vector<double>& newMemory){
  if (newMemory.size() != num_gp){
    std::ostringstream os;
    os << "Cannot set memory to new memory. Size of new Memory (" 
            << newMemory.size() << ") is not equal the number of grid points ("
            << num_gp << ") in potential with id " << id << "!";
    io::messages.add(os.str(), "BS_Potential", io::message::error);
    return;
  }
  m_auxiliaryMemory = newMemory;
}

void util::BS_Potential::getAuxMemory(std::vector<double>& newMemory){
  newMemory = m_auxiliaryMemory;
}

void 
util::BS_Potential::setMemoryParameters(double forceIncrement, int localCutoff)
{
  m_forceIncrement = forceIncrement;
  m_localCutoff = localCutoff;
}

std::string util::BS_Potential::traj_str(){
  std::ostringstream os;
  os << std::setw(3) << id << " "
          << std::setw(3) << (m_potentialType == BS_Potential::bs_sphere ? 1 : 0) << " "
          << std::setw(12) << m_potential << " "
          << std::setw(8) << m_weight;;
  
  return os.str();
}

/********************************************************
  THE SPHERE
 */
double util::BS_Sphere::calcPotential(BS_Vector & bs_pos){
  // vec(r_k)
  DEBUG(8, "\nSphere " << id << " Center: " << m_center.str() << " Position: " << bs_pos.str());
  BS_Vector radialVector; 
  bs_pos.minus(m_center, radialVector);
  double radialDistance = radialVector.normalize();
  DEBUG(6, "; RadialVector |" << radialVector.str() << "| = " << radialDistance);
  
  if (radialDistance >= m_radius){
    double diff2 = radialDistance - m_radius;
    m_potDerivatives = radialVector * (diff2 * m_force_const);
    DEBUG(10, "Multiply radial Vector with diff^2 * force = " << (diff2 * m_force_const));
    diff2 = diff2 * diff2;
    m_potential = m_memory[num_gp - 1];
    m_potential += m_half_force_const * diff2;
    m_activeGridPoint = (num_gp - 1);
    DEBUG(8, "Outside: Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Sphere " << id << " Derivatives = " << m_potDerivatives.str()); 
  }
  else {
    double gridPointCoordinate = m_gridPointScaling * radialDistance;
    //double gridPointCoordinate = (num_gp - 1) * pow((radialDistance / m_radius), radialVector.size());
    DEBUG(6, "GridPointCoordinate = " << gridPointCoordinate);
    double force;
    calcPotAndForce(gridPointCoordinate, m_potential, force);
    DEBUG(10, "Force = " << force);
    m_potDerivatives = radialVector * (m_gridPointScaling * force);
    //m_potDerivatives = radialVector * (gridPointCoordinate * radialVector.size() / radialDistance * force);
    DEBUG(6, "Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Sphere " << id << " Derivatives = " << m_potDerivatives.str());
  }
  return m_potential;
}

std::string util::BS_Sphere::str(){
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
double util::BS_Stick::calcPotential(BS_Vector& bs_pos){
  BS_Vector relativePosition;
  bs_pos.minus(m_startPoint, relativePosition);
  double longitudinalDistance = relativePosition.dot(m_unitLongitudinal);
  DEBUG(8, "\nStick " << id << "; relativePosition " << relativePosition.str() << "longDist " << longitudinalDistance);
  DEBUG(8, "\tUnit Longitudinal Vector: " << m_unitLongitudinal.str());
  // At the beginning of the stick
  if (longitudinalDistance <= 0){
    double radialDistance = relativePosition.normalize();
    m_activeGridPoint = 0;
    m_potential = m_memory[0];
    double offset = radialDistance - m_half_width;
    DEBUG(8, "offset = " << offset);
    if (offset > 0){
      m_potential += m_half_force_const * offset * offset;
      m_potDerivatives = relativePosition * (offset * m_force_const);
    } 
    else {
      m_potDerivatives.assign(relativePosition.size(), 0);
    }
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Stick " << id << " Derivatives = " << m_potDerivatives.str());
    //return false;
  }
  // At the end of the stick
  else if (longitudinalDistance >= m_length){
    bs_pos.minus(m_endPoint, relativePosition);
    double radialDistance = relativePosition.normalize();
    m_activeGridPoint = num_gp - 1;
    m_potential = m_memory[num_gp - 1];
    double offset = radialDistance - m_half_width;
    DEBUG(8, "offset = " << offset);
    if (offset > 0){
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
  else{
    double gridPointCoordinate = m_gridPointScaling * longitudinalDistance;
    double force;
    calcPotAndForce(gridPointCoordinate, m_potential, force);
    m_potDerivatives = m_unitLongitudinal * (force * m_gridPointScaling);
    BS_Vector transversalVector;
    double transversalDistance;
    relativePosition.minus(m_unitLongitudinal * longitudinalDistance,
                           transversalVector);
    transversalDistance = transversalVector.normalize();
    DEBUG(8, "Transversal distance: " << transversalDistance << "; from : " << transversalVector.str());
    double offset = transversalDistance - m_half_width;
    if (offset > 0){
      m_potential += m_half_force_const * offset * offset;
      m_potDerivatives += transversalVector * (offset * m_force_const);
    }
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << m_activeGridPoint);
    DEBUG(8, "Stick " << id << " Derivatives = " << m_potDerivatives.str());
  }   
  return m_potential;
}

std::string util::BS_Stick::str(){
  std::ostringstream os;
  os << std::setw(3) << id
          << std::setw(10) << m_force_const << " "
          << std::setw(6) << m_half_width << " "
          << m_startPoint.str() << " "
          << m_endPoint.str();
  return os.str();
}
