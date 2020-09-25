/**
 * @file shake.cc
 * contains the template methods for
 * the class Shake.
 */
#ifdef XXMPI
#include <mpi.h>
#endif

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../math/periodicity.h"

#include "../../algorithm/constraints/shake.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"
#include "../../util/undirected_graph.h"

#include <limits>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

// dihedral constraints template method
#include "dihedral_constraint.cc"

/**
 * Constructor.
 */
algorithm::Shake
::Shake(double const solute_tolerance, double const solvent_tolerance,
        int const max_iterations, std::string const name)
: Algorithm(name),
m_solute_tolerance(solute_tolerance),
m_solvent_tolerance(solvent_tolerance),
m_max_iterations(max_iterations) {
}

/**
 * Destructor.
 */
algorithm::Shake
::~Shake() {
}

void algorithm::Shake
::solute_tolerance(double const tol) {
  m_solute_tolerance = tol;
}

void algorithm::Shake
::solvent_tolerance(double const tol) {
  m_solvent_tolerance = tol;
}


/**
 * apply the SHAKE algorithm
 */
int algorithm::Shake::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(7, "applying SHAKE");
  if (!sim.mpi || m_rank == 0)
    m_timer.start();

  int error = 0;

  // set the constraint force to zero
  std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
          to = constrained_atoms().end();
  for (; it != to; ++it) {
    conf.old().constraint_force(*it) = 0.0;
  }

#ifdef XXMPI
  math::VArray & pos = conf.current().pos;
  if (sim.mpi) {
    // broadcast current and old coordinates and pos.

    MPI::COMM_WORLD.Bcast(&pos(0)(0), pos.size() * 3, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&conf.old().pos(0)(0), conf.old().pos.size() * 3, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&conf.current().box(0)(0), 9, MPI::DOUBLE, 0);

    // set virial tensor and solute coordinates of slaves to zero
    if (m_rank) { // slave
      conf.old().virial_tensor = 0.0;
      conf.old().constraint_force = 0.0;
    }
    // set positions of constraint groups worked on by other ranks to 0
    for(unsigned int i = 0; i < m_constraint_groups[m_rank].unaffected_indices.size(); ++i){
      pos(m_constraint_groups[m_rank].unaffected_indices[i]) = 0;
    }
  }
#endif

  // check whether we shake
  if ((topo.solute().distance_constraints().size() &&
          sim.param().constraint.solute.algorithm == simulation::constr_shake &&
          sim.param().constraint.ntc > 1) ||
          sim.param().dihrest.dihrest == simulation::dihedral_constr) {

    DEBUG(8, "\twe need to shake SOLUTE");

    SPLIT_VIRIAL_BOUNDARY(solute,
            topo, conf, sim,
            m_max_iterations, error);

    if (error) {
      std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLUTE "
              << "at step " << sim.steps() << std::endl;
      conf.special().shake_failure_occurred = true;
      if (!sim.mpi || m_rank == 0)
        m_timer.stop();
      return E_SHAKE_FAILURE_SOLUTE;
    }
  }

#ifdef XXMPI
  int transfer_size=topo.num_atoms();
#endif
  if (!error && sim.param().system.nsm &&
          sim.param().constraint.solvent.algorithm == simulation::constr_shake) {

    DEBUG(8, "\twe need to shake SOLVENT");
    SPLIT_VIRIAL_BOUNDARY(solvent,
            topo, conf, sim, sim.time_step_size(),
            m_max_iterations, error);
    if (error) {
      std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLVENT "
              << "at step " << sim.steps() << std::endl;
      conf.special().shake_failure_occurred = true;
      if (!sim.mpi || m_rank == 0)
        m_timer.stop();
      return E_SHAKE_FAILURE_SOLVENT;
    }
  } else {
#ifdef XXMPI
    // sum up only solute positions in the next block, otherwise things will 
    // go wrong, because all ranks still have all solvent positions
    transfer_size=topo.num_solute_atoms();
#endif
  }
#ifdef XXMPI
  if (sim.mpi) {
    if (m_rank == 0) {
      // Master 
      // reduce current positions, store them in new_pos and assign them to current positions
      math::VArray new_pos=pos;
      MPI::COMM_WORLD.Reduce(&pos(0)(0), &new_pos(0)(0),
              transfer_size * 3, MPI::DOUBLE, MPI::SUM, 0);
      pos = new_pos;

      // reduce current virial tensor, store it in virial_new and reduce it to current tensor
      math::Matrix virial_new(0.0);
      MPI::COMM_WORLD.Reduce(&conf.old().virial_tensor(0, 0), &virial_new(0, 0),
              9, MPI::DOUBLE, MPI::SUM, 0);
      conf.old().virial_tensor = virial_new;

      // reduce current contraint force, store it in cons_force_new and reduce
      // it to the current constraint force
      math::VArray cons_force_new=conf.old().constraint_force;
      MPI::COMM_WORLD.Reduce(&conf.old().constraint_force(0)(0), &cons_force_new(0)(0),
              transfer_size * 3, MPI::DOUBLE, MPI::SUM, 0);
      conf.old().constraint_force = cons_force_new;
    } else {
      // slave
      // reduce pos
      MPI::COMM_WORLD.Reduce(&pos(0)(0), NULL,
              transfer_size * 3, MPI::DOUBLE, MPI::SUM, 0);
      // reduce virial
      MPI::COMM_WORLD.Reduce(&conf.old().virial_tensor(0, 0), NULL,
              9, MPI::DOUBLE, MPI::SUM, 0);
      // reduce constraint force
      MPI::COMM_WORLD.Reduce(&conf.old().constraint_force(0)(0), NULL,
              transfer_size * 3, MPI::DOUBLE, MPI::SUM, 0);
    }
  }
#endif

  // shaken velocity:
  // stochastic dynamics, energy minimisation, analysis needs to shake without
  // velocity correction (once; it shakes twice...)
  if (!sim.param().stochastic.sd && !sim.param().minimise.ntem &&
          !sim.param().analyze.analyze) {
    std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.current().vel(*it) = (conf.current().pos(*it) - conf.old().pos(*it)) /
              sim.time_step_size();
    }
  }

  if (!sim.mpi || m_rank == 0)
    m_timer.stop();
  // return success!
  return 0;

}

int algorithm::Shake::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) {
  if (!quiet) {
    os << "SHAKE\n"
            << "\tsolute\t";
    if (sim.param().constraint.solute.algorithm == simulation::constr_shake) {
      os << "ON\n";
      os << "\t\ttolerance = "
              << sim.param().constraint.solute.shake_tolerance << "\n";
    } else os << "OFF\n";

    os << "\tsolvent\t";

    if (sim.param().constraint.solvent.algorithm == simulation::constr_shake) {
      if (sim.mpi)
        os << "ON (MPI parallel version)\n";
      else
        os << "ON\n";
      os << "\t\ttolerance = "
              << sim.param().constraint.solvent.shake_tolerance << "\n";
    } else os << "OFF\n";
  }

#ifdef XXMPI
  if (sim.mpi) {
    m_rank = MPI::COMM_WORLD.Get_rank();
    m_size = MPI::COMM_WORLD.Get_size();
  } else {
    m_rank = 0;
    m_size = 1;
  }
#else
  m_rank = 0;
  m_size = 1;
#endif

  util::UndirectedGraph g(topo.num_solute_atoms());

  if (sim.param().constraint.solute.algorithm == simulation::constr_shake) {
    // loop over the constraints to find which atoms are constrained
    {
      std::vector<topology::two_body_term_struct>::const_iterator
      it = topo.solute().distance_constraints().begin(),
              to = topo.solute().distance_constraints().end();
      for (; it != to; ++it) {
        constrained_atoms().insert(it->i);
        constrained_atoms().insert(it->j);
        g.add_edge(it->i, it->j);
      }
    }
    // also add the dihedral constrained atoms
    if (sim.param().dihrest.dihrest == simulation::dihedral_constr) {
      std::vector<topology::dihedral_restraint_struct>::const_iterator
      it = topo.dihedral_restraints().begin(),
              to = topo.dihedral_restraints().end();
      for (; it != to; ++it) {
        constrained_atoms().insert(it->i);
        constrained_atoms().insert(it->j);
        constrained_atoms().insert(it->k);
        constrained_atoms().insert(it->l);
        g.add_edge(it->i, it->j);
        g.add_edge(it->k, it->j);
        g.add_edge(it->k, it->l);
      }
    }
  }
  // detect atoms connected by constraints
  util::UndirectedGraph::component_cont_t component = g.connected_components();
  // set number of constraint groups to number of MPI processes
  m_constraint_groups.resize(m_size);

  std::vector<std::set<unsigned int> > affected_indices(m_size);

if ((topo.solute().distance_constraints().size() &&
          sim.param().constraint.solute.algorithm == simulation::constr_shake &&
          sim.param().constraint.ntc > 1)) {
  const std::vector<topology::two_body_term_struct> & dist_constr = 
  topo.solute().distance_constraints();
  // assemble distance constraints into groups
  for(std::vector<topology::two_body_term_struct>::const_iterator it = dist_constr.begin();
    it != dist_constr.end(); ++it){
	// each rank should get an equal amount of constraint groups
	const unsigned int group_id = component[it->i] % m_size;
    m_constraint_groups[group_id].distance_restraints.push_back(*it);
#ifdef XXMPI
    affected_indices[group_id].insert(it->i);
    affected_indices[group_id].insert(it->j);
#endif
  }
}
  if (sim.param().dihrest.dihrest == simulation::dihedral_constr) {
	// assemble dihedral constraints into groups
    for(std::vector<topology::dihedral_restraint_struct>::const_iterator it =
      topo.dihedral_restraints().begin(); it != topo.dihedral_restraints().end(); ++it){
	  const unsigned int group_id = component[it->i] % m_size;
      m_constraint_groups[group_id].dihedral_restraints.push_back(*it);
#ifdef XXMPI
      affected_indices[group_id].insert(it->i);
      affected_indices[group_id].insert(it->j);
      affected_indices[group_id].insert(it->k);
      affected_indices[group_id].insert(it->l);
#endif
    }
  }

#ifdef XXMPI
  // set unaffected solute indices for each group
  for (unsigned int i = 0; i < topo.num_solute_atoms(); ++i) {
	  if(affected_indices[m_rank].find(i) == affected_indices[m_rank].end())
		  m_constraint_groups[m_rank].unaffected_indices.push_back(i);
  }
  // detect unconstrained solute atoms and remove from unaffected list for master process
  if(m_rank == 0) {
      std::vector<unsigned int>::iterator it_end = m_constraint_groups[0].unaffected_indices.end();
	  for (unsigned int i = 0; i < topo.num_solute_atoms(); ++i) {
		bool is_constrained = false;
		for(int group_id = 1; group_id < m_size; ++group_id) {
			is_constrained |= affected_indices[group_id].find(i) != affected_indices[group_id].end();
		}
		if(!is_constrained) {
          it_end = std::remove(m_constraint_groups[0].unaffected_indices.begin(),
                      m_constraint_groups[0].unaffected_indices.end(), i);
        }
	  }
      m_constraint_groups[0].unaffected_indices.erase(it_end, m_constraint_groups[0].unaffected_indices.end());
  }
#endif

  if (sim.param().constraint.solvent.algorithm == simulation::constr_shake) {
    // this means that all solvent atoms are constrained. Add the to the list
    for (unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) {
      constrained_atoms().insert(i);
    }
  }


  if (sim.param().start.shake_pos) {
    if (!quiet)
      os << "\n\tshaking initial positions\n";

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
    } // if shake vel
  } else if (sim.param().start.shake_vel) {
    io::messages.add("shaking velocities without shaking positions illegal.",
            "shake", io::message::error);
  }

  if (!quiet)
    os << "END\n";

  return 0;
}


