/**
 * @file m_shake.cc
 * contains the template methods for
 * the class M_Shake.
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

#include "../../algorithm/constraints/m_shake.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"
#include <cassert>

#include <limits>
#include <vector>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

#ifdef OMP
#include <omp.h>
#endif

/**
 * Constructor.
 */
algorithm::M_Shake
::M_Shake(double const tolerance, int const max_iterations,
        std::string const name)
: Algorithm(name),
m_tolerance(tolerance),
m_max_iterations(max_iterations) {
}

/**
 * Destructor.
 */
algorithm::M_Shake
::~M_Shake() {
}

void algorithm::M_Shake
::tolerance(double const tol) {
  m_tolerance = tol;
}

/**
 * do one iteration
 */
inline
int algorithm::M_Shake::m_shake_molecule
(
        configuration::Configuration & conf,
        bool & convergence,
        int first,
        double dt2i,
        const std::vector<topology::two_body_term_struct> & constr,
        math::GenericVec<math::Vec> const & dist_old,
        bool do_virial
        ) {
  convergence = true;

  // calculate constraints: solve A*l = c by direct matrix inversion

  // to store the positions and distances
  math::VArray & new_unconstr_pos = conf.current().pos;
  math::VArray & constraint_force = conf.old().constraint_force;

  // get the new distances and calculate the difference
  const math::GenericVec<math::Vec> dist_new(
          new_unconstr_pos(first + constr[0].i) - new_unconstr_pos(first + constr[0].j),
          new_unconstr_pos(first + constr[1].i) - new_unconstr_pos(first + constr[1].j),
          new_unconstr_pos(first + constr[2].i) - new_unconstr_pos(first + constr[2].j));
  const math::Vec dist2(abs2(dist_new[0]), abs2(dist_new[1]), abs2(dist_new[2]));
  const math::Vec diff(constr_length2 - dist2);


  // do we have to shake?
  const double tol2 = tolerance() * 2.0;
  if (fabs(diff(0)) >= constr_length2(0) * tol2 ||
      fabs(diff(1)) >= constr_length2(1) * tol2 ||
      fabs(diff(2)) >= constr_length2(2) * tol2) {

    // matrix A
    // to store matrix A
    math::Matrix A;
    A(0, 0) = dot(dist_old(0), dist_new(0)) * factor(0, 0);
    A(0, 1) = dot(dist_old(1), dist_new(0)) * factor(0, 1);
    A(0, 2) = dot(dist_old(2), dist_new(0)) * factor(0, 2);
    A(1, 0) = dot(dist_old(0), dist_new(1)) * factor(1, 0);
    A(1, 1) = dot(dist_old(1), dist_new(1)) * factor(1, 1);
    A(1, 2) = dot(dist_old(2), dist_new(1)) * factor(1, 2);
    A(2, 0) = dot(dist_old(0), dist_new(2)) * factor(2, 0);
    A(2, 1) = dot(dist_old(1), dist_new(2)) * factor(2, 1);
    A(2, 2) = dot(dist_old(2), dist_new(2)) * factor(2, 2);

    // vectors orthogonal? -> SHAKE error
    if (A(0, 0) < constr_length2(0) * math::epsilon ||
            A(1, 1) < constr_length2(1) * math::epsilon ||
            A(2, 2) < constr_length2(2) * math::epsilon) {
      std::cout << "M_SHAKE ERROR\n"
              << "\tatom i    : " << first + 1 << "\n"
              << "\tatom j    : " << first + 2 << "\n"
              << "\tatom k    : " << first + 3 << "\n"
              << "\tfree i    : " << math::v2s(new_unconstr_pos(first)) << "\n"
              << "\tfree j    : " << math::v2s(new_unconstr_pos(first + 1)) << "\n"
              << "\tfree k    : " << math::v2s(new_unconstr_pos(first + 2)) << "\n"
              << "\tref r_ij     : " << math::v2s(dist_old(0)) << "\n"
              << "\tref r_ik     : " << math::v2s(dist_old(1)) << "\n"
              << "\tref r_jk     : " << math::v2s(dist_old(2)) << "\n"
              << "\tr_ij         : " << math::v2s(dist_new(0)) << "\n"
              << "\tr_ik         : " << math::v2s(dist_new(1)) << "\n"
              << "\tr_jk         : " << math::v2s(dist_new(2)) << "\n"
              << "\tA(i,i)        : " << A(0, 0) << "\n"
              << "\tA(j,j)        : " << A(1, 1) << "\n"
              << "\tA(k,k)        : " << A(2, 2) << "\n"
              << "\tconstr_ij    : " << constr_length2(0) << "\n"
              << "\tconstr_ik    : " << constr_length2(1) << "\n"
              << "\tconstr_jk    : " << constr_length2(2) << "\n"
              << "\tdiff_ij      : " << diff(0) << "\n"
              << "\tdiff_ik      : " << diff(1) << "\n"
              << "\tdiff_jk      : " << diff(2) << "\n";
      return E_SHAKE_FAILURE;
    }
    // inverse
    A = inverse(A);

    const double f0 = (A(0, 0) * diff(0) + A(0, 1) * diff(1) + A(0, 2) * diff(2)) * 0.5;
    const double f1 = (A(1, 0) * diff(0) + A(1, 1) * diff(1) + A(1, 2) * diff(2)) * 0.5;
    const double f2 = (A(2, 0) * diff(0) + A(2, 1) * diff(1) + A(2, 2) * diff(2)) * 0.5;

    DEBUG(10, "Lagrange multiplier of contraint 0: " << f0);
    DEBUG(10, "Lagrange multiplier of contraint 1: " << f1);
    DEBUG(10, "Lagrange multiplier of contraint 2: " << f2);

    const math::Vec f01 = f0 * dist_old(0) + f1 * dist_old(1);
    const math::Vec f02 = f2 * dist_old(2) - f0 * dist_old(0);
    const math::Vec f12 = f1 * dist_old(1) + f2 * dist_old(2);

    // add to the forces

    constraint_force(first) += f01;
    constraint_force(first + 1) += f02;
    constraint_force(first + 2) -= f12;

    // virial
    if (do_virial) {
      for (int a = 0; a < 3; ++a) {
        for (int aa = 0; aa < 3; ++aa) {
          conf.old().virial_tensor(a, aa) -=
                  (dist_old(0)(a) * dist_old(0)(aa) * f0 +
                  dist_old(1)(a) * dist_old(1)(aa) * f1 +
                  dist_old(2)(a) * dist_old(2)(aa) * f2) * dt2i;
        }
      }
      DEBUG(12, "\tatomic virial done");
    }

    // calculate the new constrained positions
    new_unconstr_pos(first) += f01 * mass_i[0];
    new_unconstr_pos(first + 1) += f02 * mass_i[1];
    new_unconstr_pos(first + 2) -= f12 * mass_i[2];

    convergence = false;
  } // if we have to shake

  return 0;

}

/**
 * shake solvent.
 */
void algorithm::M_Shake
::solvent(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        double dt, int const max_iterations,
        int & error) {

  DEBUG(8, "\tshaking SOLVENT - M_SHAKE");
  const double dt2 = dt * dt;
  const double dt2i = 1.0 / dt2;

  if (!sim.mpi || m_rank == 0)
    m_timer.start("solvent");

  // the first atom of a solvent
  unsigned int first = topo.num_solute_atoms();
  int tot_iterations = 0;

  error = 0;
  int my_error = error;

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
      for (unsigned int i = 0; i < first; ++i) {
        pos(i) = 0.0;
      }
    }
  }
#endif

  // for all solvents
  const unsigned int num_solvent_atoms = topo.solvent(0).num_atoms();
  const bool do_virial = sim.param().pcouple.virial == math::atomic_virial;

#ifdef OMP
#pragma omp parallel 
  {
    unsigned int strid = omp_get_num_threads();
    unsigned int tid = omp_get_thread_num();

    // the first atom of a solvent
    unsigned int first = topo.num_solute_atoms() + tid * num_solvent_atoms;
#else
  {

    unsigned int strid = 1;
    unsigned int tid = 0;
#endif
    const math::VArray & old_pos = conf.old().pos;
    math::VArray & constraint_force = conf.old().constraint_force;
    // loop over the molecules
    for (unsigned int nm = tid; nm < topo.num_solvent_molecules(0) /*- num_solvent_atoms * strid*/;
            nm += strid, first += num_solvent_atoms * strid) {
      DEBUG(15, "Thread " << tid << " entering for loop");

#ifdef XXMPI
      if (sim.mpi) {
        int stride = nm + m_rank;
        DEBUG(12, "rank: " << m_rank << " nm: " << nm << " stride: " << stride);
        if (stride % m_size != 0) {
          // set current coordinates to zero.
          for (unsigned int a = 0; a < num_solvent_atoms; ++a) {
            pos(a + first) = 0.0;
          }
          // do next molecule
          continue;
        }
      }
#endif

      const std::vector<topology::two_body_term_struct> & dist_const = topo.solvent(0).distance_constraints();
      const math::GenericVec<math::Vec> dist_old(
              old_pos(first + dist_const[0].i) - old_pos(first + dist_const[0].j),
              old_pos(first + dist_const[1].i) - old_pos(first + dist_const[1].j),
              old_pos(first + dist_const[2].i) - old_pos(first + dist_const[2].j)
              );

      int num_iterations = 0;
      bool convergence = false;
      while (!convergence) {
        DEBUG(10, "\titeration" << std::setw(10) << num_iterations
                << " max_iterations" << std::setw(10) << max_iterations);

        if (m_shake_molecule(conf, convergence, first, dt2i, dist_const, dist_old, do_virial)) {

          io::messages.add("M_SHAKE error. vectors orthogonal",
                  "M_Shake::solvent", io::message::error);

          std::cout << "M_SHAKE failure in solvent!" << std::endl;
          my_error = E_SHAKE_FAILURE_SOLVENT;
          break;
        }

        //std::cout << num_iterations+1 << std::endl;
        if (++num_iterations > max_iterations) {
          io::messages.add("M_SHAKE error: too many iterations",
                  "M_Shake::solvent",
                  io::message::critical);
          my_error = E_SHAKE_FAILURE_SOLVENT;

          break;
        }
      } // while(!convergence)
      if (my_error != error) break;

      tot_iterations += num_iterations;

      constraint_force(first) *= dt2i;
      DEBUG(5, "constraint_force " << math::v2s(constraint_force(first)));
      constraint_force(first + 1) *= dt2i;
      DEBUG(5, "constraint_force " << math::v2s(constraint_force(first + 1)));
      constraint_force(first + 2) *= dt2i;
      DEBUG(5, "constraint_force " << math::v2s(constraint_force(first + 2)));
    } // molecules
  } // OMP

  // reduce everything
#ifdef XXMPI
  if (sim.mpi) {
    if (m_rank == 0) {
      // Master
      // reduce the error to all processors
      MPI::COMM_WORLD.Allreduce(&my_error, &error, 1, MPI::INT, MPI::MAX);
      //
      // reduce current positions, store them in new_pos and assign them to current positions
      math::VArray new_pos(topo.num_atoms(), math::Vec(0.0));
      MPI::COMM_WORLD.Reduce(&pos(0)(0), &new_pos(0)(0),
              pos.size() * 3, MPI::DOUBLE, MPI::SUM, 0);
      pos = new_pos;

      // reduce current virial tensor, store it in virial_new and reduce it to current tensor
      math::Matrix virial_new(0.0);
      MPI::COMM_WORLD.Reduce(&conf.old().virial_tensor(0, 0), &virial_new(0, 0),
              9, MPI::DOUBLE, MPI::SUM, 0);
      conf.old().virial_tensor = virial_new;

      // reduce current contraint force, store it in cons_force_new and reduce
      //it to the current constraint force
      math::VArray cons_force_new(topo.num_atoms(), math::Vec(0.0));
      MPI::COMM_WORLD.Reduce(&conf.old().constraint_force(0)(0), &cons_force_new(0)(0),
              conf.old().constraint_force.size() * 3, MPI::DOUBLE, MPI::SUM, 0);
      conf.old().constraint_force = cons_force_new;
    } else {
      // reduce the error to all processors
      MPI::COMM_WORLD.Allreduce(&my_error, &error, 1, MPI::INT, MPI::MAX);

      // slave
      // reduce pos
      MPI::COMM_WORLD.Reduce(&pos(0)(0), NULL,
              pos.size() * 3, MPI::DOUBLE, MPI::SUM, 0);
      // reduce virial
      MPI::COMM_WORLD.Reduce(&conf.old().virial_tensor(0, 0), NULL,
              9, MPI::DOUBLE, MPI::SUM, 0);
      // reduce constraint force
      MPI::COMM_WORLD.Reduce(&conf.old().constraint_force(0)(0), NULL,
              conf.old().constraint_force.size() * 3, MPI::DOUBLE, MPI::SUM, 0);
    }
  }
#else
  error = my_error;
#endif
  error = my_error;

  if (!sim.mpi || m_rank == 0)
    m_timer.stop("solvent");
  DEBUG(3, "total shake solvent iterations: " << tot_iterations);
} // shake solvent

/**
 * apply the M_SHAKE algorithm
 */
int algorithm::M_Shake::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(7, "applying M_SHAKE");
  if (!sim.mpi || m_rank == 0)
    m_timer.start();

  int error = 0;

  // set the constraint force to zero
  const unsigned int num_atoms = topo.num_atoms();
  for(unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
    conf.old().constraint_force(i) = 0.0;
  }

  // check whether we shake
  if (sim.param().system.nsm &&
          sim.param().constraint.solvent.algorithm == simulation::constr_m_shake) {

    DEBUG(8, "\twe need to shake SOLVENT");
    solvent(topo, conf, sim, sim.time_step_size(),
            m_max_iterations, error);
    if (error) {
      std::cout << "M_SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLVENT "
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
    const double dti = 1.0 / sim.time_step_size();
    int num_atoms = topo.num_atoms();
#ifdef OMP
    #pragma omp parallel for
#endif
    for(int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
      conf.current().vel(i) = (conf.current().pos(i) - conf.old().pos(i)) * dti;
    }
  }

  if (!sim.mpi || m_rank == 0)
    m_timer.stop();
  // return success!
  return 0;

}

int algorithm::M_Shake::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) {
  // do some important checks
  if (sim.param().posrest.posrest == simulation::posrest_const) {
    std::vector<topology::position_restraint_struct>::const_iterator
    it = topo.position_restraints().begin(), to = topo.position_restraints().end();
    for (; it != to; ++it) {
      if (it->seq > topo.num_solute_atoms()) {
        io::messages.add("position constraints in solvent do not work with M-SHAKE",
                "m_shake", io::message::error);
        return 1;
      }
    }
  }

  if (topo.num_solvents() != 1) {
    io::messages.add("M-SHAKE is not implemented for multiple solvents.",
                "m_shake", io::message::error);
    return 1;
  }
  if (topo.solvent(0).distance_constraints().size() != 3) {
    io::messages.add("M-SHAKE is only implemented for 3 solvent constraints.",
                "m_shake", io::message::error);
    return 1;
  }
  if (topo.solvent(0).atoms().size() != 3) {
    io::messages.add("M-SHAKE is only implemented for 3 solvent atoms.",
                "m_shake", io::message::error);
    return 1;
  }

  if (!quiet) {
    os << "M-SHAKE\n"
            << "\tsolvent\t";

    if (sim.param().constraint.solvent.algorithm == simulation::constr_m_shake) {
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

  // initialize some quantities
  for(unsigned int i = 0; i < 3; ++i) {
    mass_i[i] = 1.0 / topo.mass(topo.num_solute_atoms() + i);
  }
  unsigned int k = 0;
  for (std::vector<topology::two_body_term_struct>::const_iterator
    it_k = topo.solvent(0).distance_constraints().begin(),
          to_k = topo.solvent(0).distance_constraints().end(); it_k != to_k; ++it_k, ++k) {

    assert(parameter().size() > it_k->type);
    constr_length2(k) = parameter()[it_k->type].r0 * parameter()[it_k->type].r0;
    unsigned int l = 0;
    for (std::vector<topology::two_body_term_struct>::const_iterator
      it_l = topo.solvent(0).distance_constraints().begin(),
            to_l = topo.solvent(0).distance_constraints().end(); it_l != to_l; ++it_l, ++l) {
      int d11 = 0, d12 = 0, d22 = 0, d21 = 0;
      if (it_k->i == it_l->i) {
        d11 = 1;
      } else if (it_k->i == it_l->j) {
        d12 = 1;
      }
      if (it_k->j == it_l->j) {
        d22 = 1;
      } else if (it_k->j == it_l->i) {
        d21 = 1;
      }
      factor(k, l) = (d11 - d12) * mass_i[it_k->i] + (d22 - d21) * mass_i[it_k->j];
    }
  }
  DEBUG(1, "Factor matrix:\n" << math::m2s(factor));

  if (sim.param().start.shake_pos) {
    if (!quiet)
      os << "\n\tshaking initial positions\n";

    // old and current pos and vel are the same for constrained atoms...
    const unsigned int num_atoms = topo.num_atoms();
    for(unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
      conf.old().pos(i) = conf.current().pos(i);
      conf.old().vel(i) = conf.current().vel(i);
    }

    // shake the current ones
    if (apply(topo, conf, sim))
      return E_SHAKE_FAILURE;

    for(unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
      // restore the velocities
      conf.current().vel(i) = conf.old().vel(i);
      // take a step back
      conf.old().pos(i) = conf.current().pos(i);
    }

    if (sim.param().start.shake_vel) {
      if (!quiet)
        os << "\tshaking initial velocities\n";

      for(unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
        conf.current().pos(i) = conf.old().pos(i) -
                sim.time_step_size() * conf.old().vel(i);
      }

      // shake again
      if (apply(topo, conf, sim))
        return E_SHAKE_FAILURE;

      for(unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
        // restore the positions
        conf.current().pos(i) = conf.old().pos(i);
        // velocities are in opposite direction (in time)
        conf.current().vel(i) = -1.0 * conf.current().vel(i);
        conf.old().vel(i) = conf.current().vel(i);
      }
    } // if shake vel
  } else if (sim.param().start.shake_vel) {
    io::messages.add("shaking velocities without shaking positions illegal.",
            "m_shake", io::message::error);
  }

  if (!quiet)
    os << "END\n";

  return 0;
}
