/**
 * @file position_constraints.cc
 * contains the template methods for
 * the Position_Constraints class.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../algorithm/constraints/position_constraints.h"

#include "../../util/debug.h"
#include <limits>

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

  math::VArray &pos   = conf.current().pos;
  math::VArray &vel   = conf.current().vel;
  math::VArray &force = conf.current().force;
  
  for( ; it != to; ++it){
    pos(it->seq) = conf.special().reference_positions(it->seq);
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
  math::VArray &vel   = conf.current().vel;

  std::set<unsigned int> constrained;
  for( ; it != to; ++it) {
    pos(it->seq) = conf.special().reference_positions(it->seq);
    vel(it->seq) = 0.0;
    topo.inverse_mass()(it->seq) = 0.0;
    constrained.insert(it->seq);
  }
  
  // Here, we have to check whether no atoms that are positionally
  // contrained are at the same time member of a distance constraint.
  
  // loop over distance constraints
  std::vector<topology::two_body_term_struct> keep;
  for(std::vector<topology::two_body_term_struct>::const_iterator
      dist_it = topo.solute().distance_constraints().begin(),
      dist_to = topo.solute().distance_constraints().end();
      dist_it != dist_to; ++dist_it) {
    // search for positinally contrained atoms in distance constraint
    if (constrained.count(dist_it->i) && constrained.count(dist_it->j)) {
      // remove
    } else {
      // keep
      keep.push_back(*dist_it);
    }
  }
  topo.solute().distance_constraints() = keep;

  return 0;
}
  

