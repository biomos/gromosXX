/**
 * @file flexible_constraint.cc
 * implements flexible constraint algorithm
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/nonbonded/interaction/nonbonded_interaction.h"

#include "../../interaction/forcefield/forcefield.h"

#include "../../math/periodicity.h"

#include "../../algorithm/constraints/flexible_constraint.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * Constructor.
 */
algorithm::Flexible_Constraint
::Flexible_Constraint(double const tolerance, int const max_iterations,
		      interaction::Forcefield * ff)
  : Algorithm("FlexibleShake"),
    m_tolerance(tolerance),
    m_max_iterations(max_iterations),
    m_flex_len()
{

  if (ff){
    for(size_t i=0; i<ff->size(); ++i){
      if ((*ff)[i]->name == "NonBonded"){
	// we have a nonbonded, try to cast it
	m_nonbonded = dynamic_cast<interaction::Nonbonded_Interaction *>((*ff)[i]);
      }
    }
    
    if(!ff){
      io::messages.add("Accessing the Nonbonded_Interaction failed!",
		       "Flexible_Constraint::Constructor",
		       io::message::error);
    }
  }
}

/**
 * Destructor
 */
algorithm::Flexible_Constraint
::~Flexible_Constraint()
{
}

void algorithm::Flexible_Constraint
::tolerance(double const tol)
{
  m_tolerance = tol;
}

/**
 * apply the Flexible SHAKE algorithm
 */
int algorithm::Flexible_Constraint
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "applying Flexible SHAKE");
  int error = 0;

  conf.special().flexible_constraint.flexible_ekin.assign
    (conf.special().flexible_constraint.flexible_ekin.size(), 0.0);

  // check whether we shake
  if (topo.solute().distance_constraints().size() && 
      sim.param().constraint.solute.algorithm == simulation::constr_flexshake &&
      sim.param().constraint.ntc > 1){
    DEBUG(8, "\twe need to flexible shake SOLUTE");

    calc_distance(topo, conf, sim);
    solute(topo, conf, sim, error);
  }
  
  if (error){
    std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLUTE "
	      << "at step " << sim.steps() << std::endl;
    conf.special().shake_failure_occurred = true;
    return E_SHAKE_FAILURE_SOLUTE;
  }

  // return success!
  return 0;
}


    

//================================================================================
//================================================================================
// 
// exact interation
//
//================================================================================
//================================================================================

    


//================================================================================
// flexible constraint distance
//================================================================================

void algorithm::Flexible_Constraint::calc_distance
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim
 )
{
  SPLIT_BOUNDARY(_calc_distance, topo, conf, sim);
}



//================================================================================
// undetermined forces
//================================================================================



/**
 * shake solute
 */
void algorithm::Flexible_Constraint::solute
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 int & error
 )
{
  SPLIT_VIRIAL_BOUNDARY(_solute,
			topo, conf, sim, error);
}




int algorithm::Flexible_Constraint
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       std::ostream & os,
       bool quiet)
{
  if (!quiet){
    os << "FLEXIBLESHAKE\n"
	      << "\tsolute\t";
    if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
      os << "ON\n";
      os << "\t\ttolerance = " 
		<< sim.param().constraint.solute.shake_tolerance << "\n";
      if (sim.param().constraint.solute.flexshake_readin)
	os << "\t\treading velocities along constraints from file\n";
      if (sim.param().constraint.solute.flexshake_mode == 2 ||
	  sim.param().constraint.solute.flexshake_mode == 3)
	os << "\t\tusing the exact algorithm\n";
      else
	os << "\t\tusing the approximate algorithm\n";
      if (sim.param().constraint.solute.flexshake_mode == 0 ||
	  sim.param().constraint.solute.flexshake_mode == 2)
	os << "\t\tusing potential and kinetic energy\n";
      else
	os << "\t\tusing potentialenergy only\n";
      
    }
    else os << "OFF\n";
  
    os << "\tsolvent\t";
  }
  
  if (sim.param().constraint.solvent.algorithm == simulation::constr_flexshake){
    if (!quiet)
      os << "not supported!\n";
    io::messages.add("flexible shake for solvent not implemented", "Flexible_Constraint",
		     io::message::error);
  }
  else if (!quiet) os << "OFF\n";

  if (!quiet)
    os << "END\n";

  if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake) {
    // loop over the constraints to find out which atoms are constrained
    std::vector<topology::two_body_term_struct>::const_iterator
    it = topo.solute().distance_constraints().begin(),
            to = topo.solute().distance_constraints().end();
    for (; it != to; ++it) {
      constrained_atoms().insert(it->i);
      constrained_atoms().insert(it->j);
    }
  }

  if (sim.param().start.shake_pos){
    if (!quiet)
      os << "(flexible) shaking initial positions\n";

   // old and current pos and vel are the same for constrained atoms...
    std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.old().pos(*it) = conf.current().pos(*it);
      conf.old().vel(*it) = conf.current().vel(*it);
    }

    // shake the current ones
    if (apply(topo, conf, sim))
      return E_SHAKE_FAILURE;


    it = constrained_atoms().begin();
    for (; it != to; ++it) {
      // restore the velocities
      conf.current().vel(*it) = conf.old().vel(*it);
      // take a step back
      conf.old().pos(*it) = conf.current().pos(*it);
    }
    
    conf.special().flexible_constraint.flexible_vel.assign(
    conf.special().flexible_constraint.flexible_vel.size(), 0.0);
    
    if (sim.param().start.shake_vel) {
      if (!quiet)
        os << "\tshaking initial velocities\n";

      it = constrained_atoms().begin();
      for (; it != to; ++it) {
        conf.current().pos(*it) = conf.old().pos(*it) -
                sim.time_step_size() * conf.old().vel(*it);
      }

      // shake again
      if (apply(topo, conf, sim))
        return E_SHAKE_FAILURE;

      it = constrained_atoms().begin();
      for (; it != to; ++it) {
        // restore the positions
        conf.current().pos(*it) = conf.old().pos(*it);
        // velocities are in opposite direction (in time)
        conf.current().vel(*it) = -1.0 * conf.current().vel(*it);
        conf.old().vel(*it) = conf.current().vel(*it);
      }
      
      for(unsigned int i = 0; i < conf.special().flexible_constraint.flexible_vel.size(); ++i) {
        conf.special().flexible_constraint.flexible_vel[i] *= -1.0;
      }
    } // if shake vel    
  }  else if (sim.param().start.shake_vel){
    io::messages.add("shaking velocities without shaking positions illegal.",
		     "shake", io::message::error);
  }

  os << "END\n"; 
  return 0;
}




void algorithm::Flexible_Constraint::_store_lengths
(configuration::Configuration & conf
) {
  if (conf.special().flexible_constraint.flex_len.size() < m_flex_len.size())
    conf.special().flexible_constraint.flex_len.resize(m_flex_len.size());

  for (unsigned int k = 0; k < m_flex_len.size(); ++k)
    conf.special().flexible_constraint.flex_len[k] = m_flex_len[k];

}

