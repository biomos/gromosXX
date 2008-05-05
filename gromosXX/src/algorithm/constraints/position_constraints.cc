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
       std::ostream & os,
       bool quiet)
{
  if (!quiet){
    os << "POSITION CONSTRAINTS\n"
       << "\tenabled\n"
       << "END\n";
  }
  
  // Set the particles to their position
  
  // loop over restraints and set the position
  std::vector<topology::position_restraint_struct>::const_iterator
          it = topo.position_restraints().begin(),
          to = topo.position_restraints().end();
  
  math::VArray &pos   = conf.current().pos;
  
  for( ; it != to; ++it)
    pos(it->seq) = it->pos;
  
  // Here, we have to check whether no atoms that are positionally
  // contrained are at the same time member of a distance constraint.
  
  // loop over distance constraints
  for(std::vector<topology::two_body_term_struct>::const_iterator
      dist_it = topo.solute().distance_constraints().begin(),
      dist_to = topo.solute().distance_constraints().end();
      dist_it != dist_to; ++dist_it) {
    // search for positinally contrained atoms in distance constraint
    for(std::vector<topology::position_restraint_struct>::const_iterator 
       pos_it = topo.position_restraints().begin(),
       pos_to = topo.position_restraints().end();
       pos_it != pos_to; ++pos_it) {
      
      unsigned int atom = pos_it->seq;
      if (dist_it->i == atom || dist_it->j == atom) {
        std::ostringstream msg;
        msg << "One of the atoms in bond " << dist_it->i+1 << "-" << dist_it->j+1
            << " is positionally constrained. This SHAKEing of positionally "
               "constraint atoms is not implemented.";
        io::messages.add(msg.str(), "Position_Contraints", io::message::error);
        return -1;
      }
    }
  }

  return 0;
}
  

