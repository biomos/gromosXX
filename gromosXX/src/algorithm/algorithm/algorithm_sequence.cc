/**
 * @file algorithm_sequence.cc
 * contains the inline functions for
 * algorithm_sequence.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../algorithm/constraints/remove_com_motion.h"

#include "../../algorithm/algorithm/algorithm_sequence.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

algorithm::Algorithm_Sequence::Algorithm_Sequence(bool clean)
  : std::vector<Algorithm *>(), clean(clean)
{
}

algorithm::Algorithm_Sequence::~Algorithm_Sequence()
{
  if (!clean) return;
  
  for(Algorithm_Sequence::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    delete *it;
  }
}

int algorithm::Algorithm_Sequence
::init(topology::Topology & topo, 
       configuration::Configuration &conf,
       simulation::Simulation &sim,
       std::ostream & os,
       bool quiet)
{
  DEBUG(5, "Algorithm_Sequence: init");
  // centre of mass removal
  if(!(sim.param().centreofmass.remove_trans ||
       sim.param().centreofmass.remove_rot)){

    algorithm::Remove_COM_Motion rcom;
    rcom.init(topo, conf, sim, os, quiet);
    rcom.apply(topo, conf, sim);
  }

  for(Algorithm_Sequence::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(7, "algorithm::init -> " << (*it)->name);
    if((*it)->init(topo, conf, sim, os, quiet)){
      os << "Algorithm_Sequence: error during initialisation of " 
	 << (*it)->name << "\n";
      
      io::messages.add("Error in algorithm sequence init",
		       "Algorithm_Sequence",
		       io::message::error);
    }
  }
  return 0;
}

int algorithm::Algorithm_Sequence
::run(topology::Topology & topo, 
      configuration::Configuration &conf,
      simulation::Simulation &sim)
{
  DEBUG(5, "Algorithm_Sequence: apply algorithm - START");

  for(Algorithm_Sequence::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    int ret = 0;
    DEBUG(7, "algorithm: " << (*it)->name);
    if((ret = (*it)->apply(topo, conf, sim))){
      DEBUG(1, "ERROR in algorithm_sequence::run : bailing out!");
      return ret;
    }
  }
  DEBUG(5, "Algorithm_Sequence: apply algorithm - DONE");
  return 0;
}

int algorithm::Algorithm_Sequence
::print_timing(std::ostream & os)
{
  os << "TIMING\n";
  
  for(Algorithm_Sequence::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    
    (*it)->print_timing(os);

  }
  
  os << "END\n";

  return 0;
}

algorithm::Algorithm * algorithm::Algorithm_Sequence::algorithm
(
 std::string name
 )
{
  for(Algorithm_Sequence::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    
    if ((*it)->name == name)
      return (*it);

  }

  return NULL;
}

void algorithm::Algorithm_Sequence::printSequence(){
      std::cout << " ALGORITHM SEQUENCE: \n";
      for(Algorithm_Sequence::iterator it = begin(), to = end();
      it != to;
      ++it){
          std::cout <<  "\t " << (*it)->name;
      }
      std::cout << "\n";
}