#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../interaction/interaction.h"

//#include "../math/periodicity.h"
//#include "../io/message.h"

//#include "../util/template_split.h"
#include "../util/debug.h"

#include "bs_potentials.h"
#include "bs_coordinate.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE bs_leus
    
void util::BS_Potential::calcPotAndForce(double gridPointCoordinate,
        double &potential,
        double &force){
  activeGridPoint = int (gridPointCoordinate + 0.5); // NINT(Gamma / R * r)
  int gridPointIndex1 = int(gridPointCoordinate);
  int gridPointIndex2 = int(gridPointCoordinate) + 1;
  double x = gridPointCoordinate - gridPointIndex1;
  double x2 = x * x;
  double x3 = x * x2;
  potential = memory[gridPointIndex1] * (1 - 3 * x2 + 2 * x3);
  potential += memory[gridPointIndex2] * (3 * x2 - 2 * x3);
  // potential = M[i] (1 - 3x^2 + 2x^3) + M[i+1] (1 - 3(x-1)^2 - 2(x-1)^3)
  force = (memory[gridPointIndex1] - memory[gridPointIndex2]) * 6 * (x2 - x); 
  DEBUG(10, "Grid Point Coordinate: " << gridPointCoordinate << "; active Grid Point " << activeGridPoint);
  DEBUG(10, "index1: " << gridPointIndex1 << "; index2: " << gridPointIndex2);
  DEBUG(10, "Potential: " << potential << "; force: " << force);
}  

double util::BS_Potential::calcBoltzmann(const double beta){
  m_BoltzmannFactor = exp(-beta * m_potential);
  DEBUG(10, "Boltzmann: " << m_BoltzmannFactor << " from Pot = " << m_potential << " (beta = " << beta << ")");
  return m_BoltzmannFactor;
}

double util::BS_Potential::calcWeight(double totalPartitionFct){
  m_weight = m_BoltzmannFactor / totalPartitionFct;
  // If the weight is too small, make it zero
  if (m_weight < math::epsilon) {
    m_weight = 0;
  }
  DEBUG(10, "Weight: " << m_weight << "; Boltzmann: " << m_BoltzmannFactor);
  return m_weight;
}

void util::BS_Potential::updateMemory(){
  DEBUG(10, "I update the memory by " << m_weight << " * " << memForce);
  memory[activeGridPoint] += m_weight * memForce;
}

void util::BS_Potential::setMemoryToZero(){
  memory.assign(num_gp, 0);
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
  memory = newMemory;
}

std::string util::BS_Potential::traj_str(){
  std::ostringstream os;
  os << std::setw(3) << id << " "
          << std::setw(10) << m_potential << " "
          << std::setw(10) << m_weight << " "
          << std::setw(3) << memory.size() << " ";
  
  std::vector<double>::iterator it = memory.begin(),
          to = memory.end();
  
  for (; it != to; it++){
    os << std::setw(10) << *it << " ";
  }
  return os.str();
}

/********************************************************
  THE SPHERE
 */
bool util::BS_Sphere::calcPotential(BS_Vector & bs_pos){
  // vec(r_k)
  DEBUG(8, "\nSphere " << id << " Center: " << center.str() << " Position: " << bs_pos.str());
  BS_Vector radialVector; 
  bs_pos.minus(center, radialVector);
  double radialDistance = radialVector.normalize();
  DEBUG(8, "; RadialVector |" << radialVector.str() << "| = " << radialDistance);
  
  if (radialDistance >= radius){
    double diff2 = radialDistance - radius;
    potDerivatives = radialVector * (diff2 * force_const);
    diff2 = diff2 * diff2;
    m_potential = memory[num_gp - 1];
    m_potential += half_force_const * diff2;
    activeGridPoint = (num_gp - 1);
    DEBUG(8, "Outside: Potential = " << m_potential << "; acitveGridPoint = " << activeGridPoint);
    DEBUG(8, "Sphere " << id << " Derivatives = " << potDerivatives.str()); 
    return false;
  }
  else {
    double gridPointCoordinate = gridPointScaling * radialDistance;
    double force;
    calcPotAndForce(gridPointCoordinate, m_potential, force);
    potDerivatives = radialVector * (gridPointScaling * force);
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << activeGridPoint);
    DEBUG(8, "Sphere " << id << " Derivatives = " << potDerivatives.str());
    return true;
  }
  return false;
}

std::string util::BS_Sphere::str(){
  std::ostringstream os;
  os << "Sphere " << id << "\n";
  os << "\tCenter: " << center.str() << "\n";
  return os.str();
}


/*************************************************
 THE STICK
 */
bool util::BS_Stick::calcPotential(BS_Vector& bs_pos){
  BS_Vector relativePosition;
  bs_pos.minus(startPoint, relativePosition);
  double longitudinalDistance = relativePosition.dot(unitLongitudinal);
  DEBUG(8, "\nStick " << id << "; relativePosition " << relativePosition.str() << "longDist " << longitudinalDistance);
  DEBUG(8, "\tUnit Longitudinal Vector: " << unitLongitudinal.str());
  // At the beginning of the stick
  if (longitudinalDistance <= 0){
    double radialDistance = relativePosition.normalize();
    activeGridPoint = 0;
    m_potential = memory[0];
    double offset = radialDistance - half_width;
    if (offset > 0){
      m_potential += half_force_const * offset * offset;
      potDerivatives = relativePosition * (offset * force_const);
    }
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << activeGridPoint);
    DEBUG(8, "Stick " << id << " Derivatives = " << potDerivatives.str());
    return false;
  }
  // At the end of the stick
  else if (longitudinalDistance >= length){
    bs_pos.minus(endPoint, relativePosition);
    double radialDistance = relativePosition.normalize();
    activeGridPoint = num_gp - 1;
    m_potential = memory[num_gp - 1];
    double offset = radialDistance - half_width;
    if (offset > 0){
      m_potential += half_force_const * offset * offset;
      potDerivatives = relativePosition * (offset * force_const);
    }
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << activeGridPoint);
    DEBUG(8, "Stick " << id << " Derivatives = " << potDerivatives.str());
    return false;
  } 
  // aside of the stick
  else{
    double gridPointCoordinate = gridPointScaling * longitudinalDistance;
    double force;
    calcPotAndForce(gridPointCoordinate, m_potential, force);
    potDerivatives = unitLongitudinal * (force * gridPointScaling);
    BS_Vector transversalVector;
    relativePosition.minus(unitLongitudinal * longitudinalDistance, 
            transversalVector);
    double transversalDistance = transversalVector.normalize();
    DEBUG(8, "Transversal distance: " << transversalDistance << "; from : " << transversalVector.str());
    double offset = transversalDistance - half_width;
    if (offset > 0){
      m_potential += half_force_const * offset * offset;
      potDerivatives += transversalVector * (offset * force_const);
    }
    DEBUG(8, "Potential = " << m_potential << "; acitveGridPoint = " << activeGridPoint);
    DEBUG(8, "Stick " << id << " Derivatives = " << potDerivatives.str());
    return true;
  }   
}

std::string util::BS_Stick::str(){
  std::ostringstream os;
  os << "Stick " << id << "\n";
  os << "\tStart: " << startPoint.str() << "\n";
  os << "\tEnd: " << endPoint.str() << "\n";
  return os.str();
}
