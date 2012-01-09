/**
 * @file bs_subspace.cc
 * A subspace of the B&S-LEUS scheme
 */
#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../math/periodicity.h"

#include "../util/bs_coordinate.h"
#include "../util/bs_vector.h"
#include "../util/bs_subspace.h"
#include "../util/bs_potentials.h"
#include "../util/template_split.h"
#include "../util/debug.h"


#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE bs_leus

void util::BS_Subspace::freeMemory(){
  DEBUG(8, "Free the memory of the subspace");
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
                                        pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    delete *pot_i;
  }
  std::vector<util::BS_Coordinate *>::iterator coord_i = m_definition.begin(),
                                          coord_end = m_definition.end();
  for (; coord_i != coord_end; coord_i++){
    delete *coord_i;
  }
}

double util::BS_Subspace::calculatePotential(configuration::Configuration& conf,
        simulation::Simulation &sim){
  DEBUG(8, "Subspace: Calculate the potential");
  DEBUG(8, this->debug_str());
  transformCurrent(conf);
  double totalPartitionFct = 0;
  const double beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);

  DEBUG(10, "beta = " << beta << "(k = " << math::k_Boltzmann << ")");
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
                                        pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    (*pot_i)->calcPotential(bs_pos);
    totalPartitionFct += (*pot_i)->calcBoltzmann(beta);
  }
  DEBUG(5, "Total Partition Function = " << totalPartitionFct);
  return totalPartitionFct;
}

void util::BS_Subspace::addForces(configuration::Configuration& conf, 
                                  double totalPartitionFct){
  DEBUG(8, "Subspace: add forces. (Total Partition Fct: " << totalPartitionFct << ")");
  //DEBUG(10, this->debug_str()); 
  m_force = bs_pos;
  m_force.nullify();
  BS_Vector derivatives;
  BS_Vector::iterator der_i;
  
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
                                        pot_end = potentials.end();
  std::vector<BS_Coordinate *>::iterator coord_i, coord_end = m_definition.end();
  // For every Potential
  for (; pot_i != pot_end; pot_i++){
    (*pot_i)->getDerivatives(derivatives);
    derivatives.scale((*pot_i)->calcWeight(totalPartitionFct));
    m_force += derivatives;
    der_i = derivatives.begin();
    //int dim = 0;

    // for every internal coordinate
    DEBUG(8, "BS_Subspace: Add forces to coordinates");
    for (coord_i = m_definition.begin(); coord_i != coord_end; coord_i++){
      int size = (*coord_i)->getDimension();
      BS_Vector derivativesPerCoord;
      for (int i = 0; i < size; i++){
        derivativesPerCoord.push_back(*(der_i++));
        //DEBUG(10, "Dimension Number " << dim);
        //derivativesPerCoord.push_back(derivatives[dim++]);
      }
      (*coord_i)->addForces(conf, derivativesPerCoord);
    }
  }
}

// ======= Memory ========
void util::BS_Subspace::updateMemory(){
  DEBUG(8, "Update the Memory\n");
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
                                        pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    (*pot_i)->updateMemory();
  }
}

void util::BS_Subspace::setMemory(int id, BS_Potential::potential_enum type, 
        std::vector<double> &memory)
{
  DEBUG(8, "Set the memory");
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
                                        pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->m_potentialType == type && (*pot_i)->id == id) {
      (*pot_i)->setMemory(memory);
      return;
    }
  }
  
  // We should not be here!!
  std::ostringstream msg;
  msg << "Could not set the memory of ";
  if (type == BS_Potential::bs_sphere)
    msg << "sphere";
  else
    msg << "stick";
  msg << " " << id << "!\n";
  io::messages.add(msg.str(), "BS_Subspace", io::message::error);
}

bool util::BS_Subspace::getMemory(int id, BS_Potential::potential_enum type, 
        std::vector<double>& memory) const
{
  std::vector<BS_Potential *>::const_iterator pot_i = potentials.begin(),
                                              pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->m_potentialType == type && (*pot_i)->id == id) {
      DEBUG(8, "Found Memory");
      (*pot_i)->getMemory(memory);
      return true;
    }
  }
  
  // We should not be here!!
  std::ostringstream msg;
  msg << "Could not get the memory of ";
  if (type == BS_Potential::bs_sphere)
    msg << "sphere";
  else
    msg << "stick";
  msg << " " << id << "!\n";
  io::messages.add(msg.str(), "BS_Subspace", io::message::error);
  
  return false;
}


void util::BS_Subspace::setMemoryToZero(){
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
          pot_end = potentials.end();
  
  for (; pot_i != pot_end; pot_i++){
    (*pot_i)->setMemoryToZero();
  }
}

// ===== Auxiliary Memory =====
void util::BS_Subspace::setAuxMemory(int id, BS_Potential::potential_enum type, 
        std::vector<double> &memory, int auxCounter, int redCounter)
{
  DEBUG(8, "Set the auxiliary memory");
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
                                        pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->m_potentialType == type && (*pot_i)->id == id) {
      (*pot_i)->setAuxMemory(memory, auxCounter, redCounter);
      return;
    }
  }
  
  // We should not be here!!
  std::ostringstream msg;
  msg << "Could not set the auxiliary memory of ";
  if (type == BS_Potential::bs_sphere)
    msg << "sphere";
  else
    msg << "stick";
  msg << " " << id << "!\n";
  io::messages.add(msg.str(), "BS_Subspace", io::message::error);
}

bool util::BS_Subspace::getAuxMemory(int id, BS_Potential::potential_enum type, 
        std::vector<double>& memory, int &auxCounter, int &redCounter) const
{
  std::vector<BS_Potential *>::const_iterator pot_i = potentials.begin(),
                                              pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->m_potentialType == type && (*pot_i)->id == id) {
      DEBUG(8, "Found auxiliary Memory");
      (*pot_i)->getAuxMemory(memory, auxCounter, redCounter);
      return true;
    }
  }
  
  // We should not be here!!
  std::ostringstream msg;
  msg << "Could not get the auxiliary memory of ";
  if (type == BS_Potential::bs_sphere)
    msg << "sphere";
  else
    msg << "stick";
  msg << " " << id << "!\n";
  io::messages.add(msg.str(), "BS_Subspace", io::message::error);
  
  return false;
}

void util::BS_Subspace::setAuxMemoryToZero(){
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
          pot_end = potentials.end();
  
  for (; pot_i != pot_end; pot_i++){
    (*pot_i)->setAuxMemoryToZero();
  }
}


// ===== Coordinates
void util::BS_Subspace::transformCurrent(configuration::Configuration& conf){
  transformCoordinates(conf, bs_pos);
}

void util::BS_Subspace::transformCoordinates(configuration::Configuration& conf, 
        BS_Vector& result){
  DEBUG(8, "Subspace: Transform the coordinates");
  result.clear();
  std::vector<BS_Coordinate *>::iterator def_i = m_definition.begin(),
          def_end = m_definition.end();
  BS_Vector coord;
  BS_Vector::iterator it, to;
  
  for (; def_i != def_end; def_i++){
    (*def_i)->calculateInternalCoord(conf);
    (*def_i)->getInternalCoordinates(coord);
    for(it = coord.begin(), to = coord.end(); it != to; it++){
      result.push_back(*it);
    }
  }
}

void util::BS_Subspace::addCoordinate(BS_Coordinate* newCoordinate){
  DEBUG(8, "Subspace: Added new coordinate")
  m_definition.push_back(newCoordinate);
}

void util::BS_Subspace::addPotential(BS_Potential* newPotential){
  newPotential->setMemoryParameters(m_forceIncrement, m_reductionFactor, 
                                    m_localCutoff, m_globalCutoff);
  potentials.push_back(newPotential);
  if (newPotential->m_potentialType == BS_Potential::bs_sphere){
    DEBUG(8, "Subspace: Added new Sphere");
    m_numSpheres++;
  }
  else if (newPotential->m_potentialType == BS_Potential::bs_stick){
    DEBUG(8, "Subspace: Added new Stick");
    m_numSticks++;
  }
  else {
    io::messages.add("Added unknown type of potential to Subspace",
            "BS_Subspace", io::message::error);
  }
}

util::BS_Vector util::BS_Subspace::getCenter(int id){
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
                                        pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->m_potentialType == BS_Potential::bs_sphere &&
            (*pot_i)->id == id) {
      return (static_cast<util::BS_Sphere *> (*pot_i))->getCenter();
    }
  }
  std::ostringstream msg;
  msg << "Could not obtain the center of sphere " << id << "!\n";
  io::messages.add(msg.str(), "BS_Subspace", io::message::error);
}


void util::BS_Subspace::getForce(std::vector<double>& force){
  if (m_force.size() == 0){
    force.assign(m_definition.size(), 0);
    return;
  }
  BS_Vector::iterator force_i = m_force.begin(),
          force_end = m_force.end();
  force.clear();
  for (; force_i != force_end; force_i++){
    force.push_back(*force_i);
  }
}


std::string util::BS_Subspace::debug_str(){
  std::ostringstream os;
  os << "Subspace\n";
  std::vector<BS_Potential *>::iterator pot_i = potentials.begin(),
                                        pot_end = potentials.end();
  for (; pot_i != pot_end; pot_i++){
    os << (*pot_i)->str();
  }
  std::vector<util::BS_Coordinate *>::iterator coord_i = m_definition.begin(),
                                          coord_end = m_definition.end();
  for (; coord_i != coord_end; coord_i++){
    os << (*coord_i)->str();
  }
  os << "Current Position: " << bs_pos.str();
  os << "\n";
  return os.str();
}

std::string util::BS_Subspace::traj_str(){
  std::ostringstream os;
  os << "# NumSpheres NumSticks\n"
          << std::setw(8) << m_numSpheres << std::setw(8) << m_numSticks << "\n";

  if (m_numSpheres > 0) {
    os << "# Spheres\n" << "# ID Potential  Weight NumMemPoints [Memory...]\n";
    std::vector<BS_Potential *>::iterator it = potentials.begin(),
            to = potentials.end();
    for (; it != to; it++) {
      if ((*it)->m_potentialType == BS_Potential::bs_sphere) {
        os << (*it)->traj_str() << "\n";
      }
    }
  }

  if (m_numSticks > 0) {
    os << "# Sticks\n" << "# ID Potential  Weight NumMemPoints [Memory...]\n";
    std::vector<BS_Potential *>::iterator it = potentials.begin(),
            to = potentials.end();
    for (; it != to; it++) {
      if ((*it)->m_potentialType == BS_Potential::bs_stick) {
        os << (*it)->traj_str() << "\n";
      }
    }
  }
  
  return os.str();
}
