/**
 * @file position_constraints.cc
 * contains the template methods for
 * the Position_Constraints class.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/constraints/position_constraints.h>

#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * Constructor.
 */
algorithm::Position_Constraints
::Position_Constraints(std::string const name)
  : Algorithm(name)
{
}

/**
 * Destructor.
 */
algorithm::Position_Constraints
::~Position_Constraints()
{
}

/**
 * apply roto-translational constraints
 */
int algorithm::Position_Constraints
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "position constraints");

  // loop over the position restraints
  std::vector<topology::position_restraint_struct>::const_iterator 
    it = topo.position_restraints().begin(),
    to = topo.position_restraints().end();

  math::VArray &vel   = conf.current().vel;
  math::VArray &force = conf.current().force;
  
  for( ; it != to; ++it){

    force(it->seq) = 0.0;
    vel(it->seq) = 0.0;

  }
  
  return 0;		   
}

int algorithm::Position_Constraints
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       bool quiet)
{
  if (!quiet){
    std::cout << "POSITION CONSTRAINTS\n"
	      << "\tenabled\n"
	      << "\tSHAKEing of bonds containing positionally constrained"
	      << " atoms is NOT implemented (results will be wrong)!\n"
	      << "END\n";
  }

  return 0;
}
  

