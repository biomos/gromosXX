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
::Shake(double const tolerance, int const max_iterations,
        std::string const name)
: Algorithm(name),
m_tolerance(tolerance),
m_max_iterations(max_iterations) {
}

/**
 * Destructor.
 */
algorithm::Shake
::~Shake() {
}

void algorithm::Shake
::tolerance(double const tol) {
  m_tolerance = tol;
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

  // check whether we shake
  if (m_rank == 0 &&
          ((topo.solute().distance_constraints().size() &&
          sim.param().constraint.solute.algorithm == simulation::constr_shake &&
          sim.param().constraint.ntc > 1) ||
          sim.param().dihrest.dihrest == simulation::dihedral_constr)) {

    DEBUG(8, "\twe need to shake SOLUTE");

    SPLIT_VIRIAL_BOUNDARY(solute,
            topo, conf, sim,
            m_max_iterations, error);

#ifdef XXMPI
    if (sim.mpi && sim.param().system.nsm &&
      sim.param().constraint.solvent.algorithm == simulation::constr_shake)
         MPI::COMM_WORLD.Bcast(&error, 1, MPI::INT, 0);
#endif
    if (error) {
      std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLUTE "
              << "at step " << sim.steps() << std::endl;
      conf.special().shake_failure_occurred = true;
      m_timer.stop();
      return E_SHAKE_FAILURE_SOLUTE;
    }
  }
#ifdef XXMPI
  //receive solute shake error
  if (sim.mpi && m_rank) {
      MPI::COMM_WORLD.Bcast(&error, 1, MPI::INT, 0);
      if (error) return E_SHAKE_FAILURE_SOLUTE;
  }
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
      m_timer.stop();
      return E_SHAKE_FAILURE_SOLVENT;
    }
  }

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

  if (sim.param().constraint.solute.algorithm == simulation::constr_shake) {
    // loop over the constraints to find out which atoms are constrained
    {
      std::vector<topology::two_body_term_struct>::const_iterator
      it = topo.solute().distance_constraints().begin(),
              to = topo.solute().distance_constraints().end();
      for (; it != to; ++it) {
        constrained_atoms().insert(it->i);
        constrained_atoms().insert(it->j);
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
      }
    }
  }

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


