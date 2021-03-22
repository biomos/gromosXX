#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../math/periodicity.h"

#include "../util/bs_umbrella.h"
#include "../util/bs_subspace.h"
#include "../util/bs_potentials.h"
#include "../util/template_split.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE bs_leus

util::BS_Umbrella::~BS_Umbrella(){
  std::vector<BS_Subspace *>::iterator it = m_subspaces.begin(),
          to = m_subspaces.end();
  for (; it !=  to; it++){
    (*it)->freeMemory();
    delete *it;
  }
}

void util::BS_Umbrella::apply(configuration::Configuration& conf,
        simulation::Simulation &sim){
  DEBUG(5, "BS_Umbrella: Apply BSLEUS");
  std::vector<double> potentials;
  std::vector<BS_Subspace *>::iterator it = m_subspaces.begin(),
          to = m_subspaces.end();
  for (; it !=  to; it++){
    DEBUG(5, "BS_Umbrella: Calculate Potentials from Subspace");
    (*it)->calculatePotential(conf, potentials);
  }
  
  double pmin = potentials[0];
  for (unsigned int i = 1; i < potentials.size(); i++){
    if (potentials[i] < pmin)
      pmin = potentials[i];
  }
  const double beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
  DEBUG(5, "BS_Umbrella: beta: " << beta);
  double sum = 0;
  for (unsigned int i = 0; i < potentials.size(); i++){
    sum += exp(-beta * (potentials[i] - pmin));
    DEBUG(8, "Potential[" << i << "] = " << potentials[i]);
  }
  m_bsleus_total = - 1.0 / beta * log(sum) + pmin;
  DEBUG(6, "Total Potential = " << m_bsleus_total);
          
  conf.current().energies.bsleus_total = m_bsleus_total;
  
  for (it = m_subspaces.begin(); it != to; it++){
    DEBUG(5, "BS_Umbrella: Add the forces from the subspace");
    (*it)->addForces(conf, potentials, beta);
  }
  
  if(sim.param().bsleus.building){
  DEBUG(5, "BS_Umbrella: Building, therefore: update the memory");
    for (it = m_subspaces.begin(); it != to; it++){
     (*it)->updateMemory();
    }
  }
}

//------------------------------------------
void util::BS_Umbrella::addSubspaces(std::vector<BS_Subspace  *> &subspaces){
  DEBUG(5, "Added the subspaces");
  m_subspaces = subspaces;
}

void util::BS_Umbrella::setCounter(int subid, int auxCounter, int redCounter)
{
  m_subspaces[subid]->setCounter(auxCounter, redCounter);
}

void util::BS_Umbrella::getCounter(int subid, 
                                   unsigned int& auxCounter, 
                                   unsigned int& redCounter) const
{
  m_subspaces[subid]->getCounter(auxCounter, redCounter);
}

// Memories
void util::BS_Umbrella::setMemory(int id, int subid, 
                                  std::vector<double>& memory)
{
  m_subspaces[subid]->setMemory(id, memory);
}

bool util::BS_Umbrella::getMemory(int id, unsigned int &subid, 
                                  std::vector<double>& memory) const
{
  for (unsigned int i = 0; i < m_subspaces.size(); i++){
    if (m_subspaces[i]->getMemory(id, memory)){
      subid = i;
      return true;
    }
  }
  return false;
}

void util::BS_Umbrella::setMemoryToZero(){
  for (unsigned int i = 0; i < m_subspaces.size(); i++)
    m_subspaces[i]->setMemoryToZero();
}

// Auxiliary Memories
void util::BS_Umbrella::setAuxMemory(int id, int subid, 
                                     std::vector<double>& memory)
{
  m_subspaces[subid]->setAuxMemory(id, memory);
}

bool util::BS_Umbrella::getAuxMemory(int id, unsigned int &subid, 
                                     std::vector<double>& memory) const
{
  for (unsigned int i = 0; i < m_subspaces.size(); i++){
    if (m_subspaces[i]->getAuxMemory(id, memory)){
      subid = i;
      return true;
    }  
  }
  return false;
}

void util::BS_Umbrella::setAuxMemoryToZero(){
  for (unsigned int i = 0; i < m_subspaces.size(); i++)
    m_subspaces[i]->setAuxMemoryToZero();
}


//------------------------------------------
void util::BS_Umbrella::getPosition(unsigned int subid, std::vector<double> &position) const
{
  assert(subid < m_subspaces.size());
  m_subspaces[subid]->getPosition(position);
}

void util::BS_Umbrella::setPosition(unsigned int subid, std::vector<double> &position)
{
  assert(subid < m_subspaces.size());
  m_subspaces[subid]->setPosition(position);
}

//------------------------------------------
double util::BS_Umbrella::getTotalPotential() const {
  return m_bsleus_total;
}

void util::BS_Umbrella::getNumPotentials(int& numPotentials) const
{
  assert(m_subspaces.size() > 0);
  numPotentials = m_subspaces[0]->getNumPotentials(); 
}

void util::BS_Umbrella::getForce(std::vector<double>& force) const {
  m_subspaces[0]->getForce(force);
}

bool util::BS_Umbrella::printAuxMem() const {
  bool printIt = true;
  for (unsigned int i = 0; i < m_subspaces.size(); i++)
    printIt = m_subspaces[i]->printAuxMem() && printIt;
  return printIt;
}

//------------------------------------------
std::string util::BS_Umbrella::traj_str() const{
  std::ostringstream os;
  os << "# BS_UMBRELLA\n"
          << "# TOT_ENERGY NUM_SUBSPACES\n"
          << std::setw(10) << m_bsleus_total << " "
          << std::setw(10) << m_subspaces.size() << "\n";
  
  os << "# BS_SUBSPACES\n";
  
  std::vector<BS_Subspace*>::const_iterator it = m_subspaces.begin(),
          to = m_subspaces.end();
  
  int i = 1;
  for (; it != to; it++, i++){
    os << "# BS_SUBSPACE " << i << "\n"
       << "# SUBID TOTDIM NUMPOTS AUXCOUNT REDCOUNT\n";
    os << (*it)->traj_str();
  }
  return os.str();
}

std::string util::BS_Umbrella::str(){
  std::ostringstream os;
  os << "# BS_UMBRELLA\n"
          << "# NUM_SUBSPACES SMOOTHING\n"
          << std::setw(10) << m_subspaces.size() << " "
          << std::setw(10) << m_smoothing << "\n";
  
  os << "# BS_SUBSPACES\n";
  
  std::vector<BS_Subspace*>::iterator it = m_subspaces.begin(),
          to = m_subspaces.end();
  
  int i = 1;
  for (; it != to; it++, i++){
    os << "# BS_SUBSPACE " << i << "\n"
       << "# SUBID NUMCOORD NUMSPH NUMSTK NUMSNK NUMPPS FORCEINCR FORCERED LOCCUTOFF GLOBCUTOFF\n";
    os << (*it)->str();
  }
  return os.str();
}
