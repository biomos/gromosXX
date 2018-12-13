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
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
                                        pot_end = m_potentials.end();
  for (; pot_i != pot_end; pot_i++){
    delete *pot_i;
  }
  std::vector<util::BS_Coordinate *>::iterator coord_i = m_definition.begin(),
                                          coord_end = m_definition.end();
  for (; coord_i != coord_end; coord_i++){
    delete *coord_i;
  }
}

void util::BS_Subspace::calculatePotential(configuration::Configuration& conf,
        std::vector<double> &potentials){
  DEBUG(8, "Subspace: Calculate the potential");
  DEBUG(8, this->debug_str());
  transformCurrent(conf);
  
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
                                        pot_end = m_potentials.end();
  for (; pot_i != pot_end; pot_i++){
    potentials.push_back((*pot_i)->calcPotential(bs_pos));
  }
}

void util::BS_Subspace::addForces(configuration::Configuration& conf, 
                                  std::vector<double> &potentials,
                                  double beta){
  m_force = bs_pos;
  m_force.nullify();
  BS_Vector derivatives;
  BS_Vector::iterator der_i;
  
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
                                        pot_end = m_potentials.end();
  std::vector<BS_Coordinate *>::iterator coord_i, coord_end = m_definition.end();
  // For every Potential
  for (; pot_i != pot_end; pot_i++){
    (*pot_i)->getDerivatives(derivatives);
    derivatives.scale((*pot_i)->calcWeight(potentials, beta));
    m_force += derivatives;
    der_i = derivatives.begin();

    // for every internal coordinate
    DEBUG(8, "BS_Subspace: Add forces to coordinates");
    for (coord_i = m_definition.begin(); coord_i != coord_end; coord_i++){
      int size = (*coord_i)->getDimension();
      BS_Vector derivativesPerCoord;
      for (int i = 0; i < size; i++){
        derivativesPerCoord.push_back(*(der_i++));
      }
      (*coord_i)->addForces(conf, derivativesPerCoord);
    }
  }
}

// ======= Memory ========
void util::BS_Subspace::updateMemory(){
  DEBUG(8, "Update the Memory");
  if (m_reductionFactor == 1.0){
    DEBUG(8, "Don't use the force reduction scheme")
    std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
                                        pot_end = m_potentials.end();
    for (; pot_i != pot_end; pot_i++) {
      (*pot_i)->updateMemory(false);
    }
  }
  else {
    DEBUG(8, "Use the force reduction scheme")
    // Are all auxilary points bigger than the local cutoff?
    bool bigger = true;
    std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
            pot_end = m_potentials.end();
    for (; pot_i != pot_end; pot_i++) {
      bigger = (*pot_i)->updateMemory(true) && bigger;
    }

    if (bigger) { // A(m,i) > local cutoff   for all m,i
      DEBUG(8, "All Auxiliary Memories were bigger than local cutoff");
      m_auxilliaryCounter++;
      bool reduceForce = false;
      if (m_auxilliaryCounter >= m_globalCutoff) {
        m_reductionCounter++;
        m_forceIncrement *= m_reductionFactor;
        m_auxilliaryCounter = 0;
        reduceForce = true;
        DEBUG(8, "Auxiliary Counter is larger than local cutoff!");
      }
      for (pot_i = m_potentials.begin(); pot_i != pot_end; pot_i++) {
        (*pot_i)->setAuxMemoryToZero();
        if (reduceForce)
            (*pot_i)->setMemoryParameters(m_forceIncrement, m_localCutoff);
      }
    }
  }
}

bool util::BS_Subspace::setMemory(int id, std::vector<double> &memory)
{
  DEBUG(8, "Set the memory");
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
                                        pot_end = m_potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->id == id) {
      (*pot_i)->setMemory(memory);
      return true;
    }
  }
  return false;
}

bool util::BS_Subspace::getMemory(int id, std::vector<double>& memory) const
{
  std::vector<BS_Potential *>::const_iterator pot_i = m_potentials.begin(),
                                              pot_end = m_potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->id == id) {
      DEBUG(8, "Found Memory");
      (*pot_i)->getMemory(memory);
      return true;
    }
  }
  return false;
}


void util::BS_Subspace::setMemoryToZero(){
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
          pot_end = m_potentials.end();
  
  for (; pot_i != pot_end; pot_i++){
    (*pot_i)->setMemoryToZero();
  }
}

// ===== Auxiliary Memory =====
bool util::BS_Subspace::setAuxMemory(int id, std::vector<double> &memory)
{
  DEBUG(8, "Set the auxiliary memory");
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
                                        pot_end = m_potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->id == id) {
      (*pot_i)->setAuxMemory(memory);
      (*pot_i)->setMemoryParameters(m_forceIncrement, m_localCutoff);
      return true;
    }
  }
  return false;
}

bool util::BS_Subspace::getAuxMemory(int id, std::vector<double>& memory) const
{
  std::vector<BS_Potential *>::const_iterator pot_i = m_potentials.begin(),
                                              pot_end = m_potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->id == id) {
      DEBUG(8, "Found auxiliary Memory");
      (*pot_i)->getAuxMemory(memory);
      return true;
    }
  }
  return false;
}

void util::BS_Subspace::setAuxMemoryToZero(){
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
          pot_end = m_potentials.end();
  
  for (; pot_i != pot_end; pot_i++){
    (*pot_i)->setAuxMemoryToZero();
  }
}


// ===== Coordinates
void util::BS_Subspace::setPosition(std::vector<double> &position){
  if (position.size() != getNumDimensions()){
    DEBUG(6, "position.size = " << position.size() << " number of Dimensions = " <<  getNumDimensions());
    io::messages.add("BS_Subspace::setPosition: position has not the right size!",
            "BS_Subspace", io::message::error);
  }
  bs_pos = position;
  std::vector<BS_Coordinate *>::iterator def_i   = m_definition.begin(),
                                         def_end = m_definition.end();
  BS_Vector::iterator it = bs_pos.begin(), to;
  
  for (; def_i != def_end; def_i++){
    to = it + (*def_i)->getDimension();
    (*def_i)->setOldPos(it, to);
    it = to;
  }
}

void util::BS_Subspace::getPosition(std::vector<double> &position){
  position = bs_pos;
}

void util::BS_Subspace::transformCurrent(configuration::Configuration& conf){
  transformCoordinates(conf, bs_pos);
  DEBUG(8, "New Position: " << bs_pos.str());
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
  newPotential->setMemoryParameters(m_forceIncrement, m_localCutoff);
  for (unsigned int i = 0; i < m_potentials.size(); i++){
    if (newPotential->id == m_potentials[i]->id){
      std::stringstream msg;
      msg << "Two Potentials of type " 
              << BS_Potential::potentialType(newPotential->m_potentialType) << " and " 
              << BS_Potential::potentialType(m_potentials[i]->m_potentialType)
              << " have the same ID (" << newPotential->id << ")!";
      io::messages.add(msg.str(), "BS_Subspace", io::message::error);
    }
  }
  m_potentials.push_back(newPotential);
  if (newPotential->m_potentialType == BS_Potential::bs_sphere){
    //DEBUG(8, "Subspace: Added new Sphere");
    m_numSpheres++;
  }
  else if (newPotential->m_potentialType == BS_Potential::bs_stick){
    //DEBUG(8, "Subspace: Added new Stick");
    m_numSticks++;
  }
  else if (newPotential->m_potentialType == BS_Potential::bs_snake){
    //DEBUG(8, "Subspace: Added new Snake");
    m_numSnakes++;
  }
  else if (newPotential->m_potentialType == BS_Potential::bs_pipe){
    //DEBUG(8, "Subspace: Added new Snake");
    m_numPipes++;
  }
  else {
    io::messages.add("Added unknown type of potential to Subspace",
            "BS_Subspace", io::message::error);
  }
  DEBUG(8, "Subspace: Added new " 
          << BS_Potential::potentialType(newPotential->m_potentialType));
}

util::BS_Vector util::BS_Subspace::getCenter(int id){
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
                                        pot_end = m_potentials.end();
  for (; pot_i != pot_end; pot_i++){
    if ((*pot_i)->m_potentialType == BS_Potential::bs_sphere &&
            (*pot_i)->id == id) {
      return (static_cast<util::BS_Sphere *> (*pot_i))->getCenter();
    }
  }
  std::ostringstream msg;
  msg << "Could not obtain the center of potential " << id << "!\n";
  io::messages.add(msg.str(), "BS_Subspace", io::message::error);
  return BS_Vector();
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

std::vector<int> util::BS_Subspace::getDimensionality() {
  std::vector<int> dimensionality;
  std::vector<BS_Coordinate *>::iterator def_i = m_definition.begin(),
                                         def_end = m_definition.end();
  
  for (; def_i != def_end; def_i++){
    dimensionality.push_back((*def_i)->getDimension());
  }
  return dimensionality;
}

unsigned int util::BS_Subspace::getNumDimensions() {
  unsigned int dimensionality = 0;
  std::vector<BS_Coordinate *>::iterator def_i = m_definition.begin(),
                                         def_end = m_definition.end();
  
  for (; def_i != def_end; def_i++){
    dimensionality += (*def_i)->getDimension();
  }
  return dimensionality;
}

bool
util::BS_Subspace::testType(int id, BS_Coordinate::Coord_type type)
{
  std::vector<BS_Coordinate *>::iterator def_i = m_definition.begin(),
                                         def_end = m_definition.end();
  
  for (; def_i != def_end; def_i++){
    if ((*def_i)->cid() == id)
      return ((*def_i)->getType() == type);
  }
  return false;
}

void
util::BS_Subspace::getCounter(unsigned int &auxilliaryCounter, 
                               unsigned int &reductionCounter)
{
  auxilliaryCounter = m_auxilliaryCounter;
  reductionCounter = m_reductionCounter;
}

void
util::BS_Subspace::setCounter(unsigned int auxilliaryCounter, 
                               unsigned int reductionCounter)
{
  m_auxilliaryCounter = auxilliaryCounter;
  m_reductionCounter = reductionCounter;
  m_forceIncrement *= pow(m_reductionFactor, m_reductionCounter);
}

// ===================================================
std::string util::BS_Subspace::debug_str(){
  std::ostringstream os;
  os << "Subspace\n";
  std::vector<BS_Potential *>::iterator pot_i = m_potentials.begin(),
                                        pot_end = m_potentials.end();
  for (; pot_i != pot_end; pot_i++){
    os << (*pot_i)->str() << "\n";
  }
  os << "\n";
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
  
  os << std::setw(3) << id << " "
     << this->getNumDimensions() << " "
     << m_potentials.size() << " "
     << m_auxilliaryCounter << " "
     << m_reductionCounter << "\n";
  
  
  os << "# Force\n"
     << m_force.str() << "\n"
     << "# Position\n"
     << bs_pos.str() << "\n";
   
  os << "# Potentials\n";
  os << "# ID TYPE   POTENTIAL     WEIGHT  ACT_GRID_POINT\n";
  std::vector<BS_Potential *>::iterator it = m_potentials.begin(),
            to = m_potentials.end();
  for (; it != to; it++) {
    os << (*it)->traj_str() << "\n";
  } 
  
  return os.str();
}

std::string
util::BS_Subspace::str()
{
  std::ostringstream os;
  os << std::setw(3) << id << " "
     << std::setw(6) << m_definition.size() << " "
     << std::setw(9) << m_numSpheres << " "
     << std::setw(7) << m_numSticks << " "
     << std::setw(7) << m_numSnakes << " "
     << std::setw(7) << m_numPipes << " "
     << std::setw(7) << m_forceIncrement << " "
     << std::setw(10) << m_reductionFactor << " "
     << std::setw(9) << m_localCutoff << " "
     << std::setw(10) << m_globalCutoff << "\n";
  
  os << "# BS_COORDINATES\n"
     << "# ID RED_FAC   DIM    TYPE\n";
  
  std::vector<BS_Coordinate*>::iterator it = m_definition.begin(),
          to = m_definition.end();
  for (; it != to; it++){
    os << (*it)->init_str() << "\n"; 
  }
  
  if (m_numSpheres) {
    os << "# BS_SPH\n"
            << "# ID FORCECST RADIUS CENTER\n";

    std::vector<BS_Potential*>::iterator pot_i = m_potentials.begin(),
            pot_end = m_potentials.end();
    for (; pot_i != pot_end; pot_i++) {
      if ((*pot_i)->m_potentialType == BS_Potential::bs_sphere)
        os << (*pot_i)->str() << "\n";
    }
  }

  if (m_numSticks) {
    os << "# BS_STK\n"
            << "# ID FORCECST WIDTH START END\n";

    std::vector<BS_Potential*>::iterator pot_i = m_potentials.begin(), 
            pot_end = m_potentials.end();
    for (; pot_i != pot_end; pot_i++) {
      if ((*pot_i)->m_potentialType == BS_Potential::bs_stick)
        os << (*pot_i)->str() << "\n";
    }
  }
  
  if (m_numSnakes) {
    os << "# BS_SNAKE\n";

    std::vector<BS_Potential*>::iterator pot_i = m_potentials.begin(), 
            pot_end = m_potentials.end();
    for (; pot_i != pot_end; pot_i++) {
      if ((*pot_i)->m_potentialType == BS_Potential::bs_snake)
        os << (*pot_i)->str() << "\n";
    }
  }
  
  if (m_numPipes) {
    os << "# BS_Pipe\n";
    os << "# ID     CLE IHW_s  OHW_s   IHW_e  OHW_e   START   END\n";

    std::vector<BS_Potential*>::iterator pot_i = m_potentials.begin(), 
            pot_end = m_potentials.end();
    for (; pot_i != pot_end; pot_i++) {
      if ((*pot_i)->m_potentialType == BS_Potential::bs_pipe)
        os << (*pot_i)->str() << "\n";
    }
  }
  return os.str();
  
}
