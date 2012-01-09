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
  double totalPartitionFct = 0; // sum[ exp (-beta *b_m) ]
  std::vector<BS_Subspace *>::iterator it = m_subspaces.begin(),
          to = m_subspaces.end();
  for (; it !=  to; it++){
    DEBUG(5, "BS_Umbrella: Calculate Potentials from Subspace");
    totalPartitionFct += (*it)->calculatePotential(conf, sim);
  }
  
  for (it = m_subspaces.begin(); it != to; it++){
    DEBUG(5, "BS_Umbrella: Add the forces from the subspace");
    (*it)->addForces(conf, totalPartitionFct);
  }
  
  const double beta_inv = (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
  DEBUG(5, "BS_Umbrella: Inverse beta: " << beta_inv)
  m_bsleus_total = - beta_inv * log(totalPartitionFct);
  conf.current().energies.bsleus_total = m_bsleus_total;
  
  if(sim.param().bsleus.building){
    for (it = m_subspaces.begin(); it != to; it++){
     DEBUG(5, "BS_Umbrella: Building, therefore: update the memory");
     (*it)->updateMemory();
    }
  }
}

//------------------------------------------
void util::BS_Umbrella::addSubspace(std::vector<BS_Subspace  *> &subspaces){
  DEBUG(5, "Added the subspaces");
  m_subspaces = subspaces;
}

// Memories
void util::BS_Umbrella::setMemory(int id, BS_Potential::potential_enum type, 
                                  std::vector<double>& memory)
{
  m_subspaces[0]->setMemory(id, type, memory);
}

bool util::BS_Umbrella::getMemory(int id, BS_Potential::potential_enum type, 
                                  std::vector<double>& memory) const
{
  return m_subspaces[0]->getMemory(id, type, memory);  
}

void util::BS_Umbrella::setMemoryToZero(){
  m_subspaces[0]->setMemoryToZero();
}

// Auxiliary Memories
void util::BS_Umbrella::setAuxMemory(int id, BS_Potential::potential_enum type, 
                                     std::vector<double>& memory, int auxCounter, 
                                     int redCounter)
{
  m_subspaces[0]->setAuxMemory(id, type, memory, auxCounter, redCounter);
}

bool util::BS_Umbrella::getAuxMemory(int id, BS_Potential::potential_enum type, 
                                     std::vector<double>& memory, int &auxCounter, 
                                     int &redCounter) const
{
  return m_subspaces[0]->getAuxMemory(id, type, memory, auxCounter, redCounter);  
}

void util::BS_Umbrella::setAuxMemoryToZero(){
  m_subspaces[0]->setAuxMemoryToZero();
}


//------------------------------------------
double util::BS_Umbrella::getTotalPotential() const {
  return m_bsleus_total;
}

void util::BS_Umbrella::getNumPotentials(int& numSpheres, int& numSticks) const
{
  numSticks = m_subspaces[0]->getNumSticks(); 
  numSpheres = m_subspaces[0]->getNumSpheres();
}

void util::BS_Umbrella::getForce(std::vector<double>& force) const {
  m_subspaces[0]->getForce(force);
}

//------------------------------------------
std::string util::BS_Umbrella::traj_str() const{
  std::ostringstream os;
  os << m_subspaces[0]->traj_str();
  return os.str();
}
