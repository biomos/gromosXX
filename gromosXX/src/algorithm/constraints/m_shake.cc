/**
 * @file m_shake.cc
 * contains the template methods for
 * the class M_Shake.
 */
/*#ifdef XXMPI
#include <mpi.h>
#endif*/

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>

#include <math/periodicity.h>

#include <algorithm/constraints/m_shake.h>

#include <util/template_split.h>
#include <util/error.h>
#include <util/debug.h>
#include <cassert>

#include <limits>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

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
template<math::boundary_enum B, math::virial_enum V>
int algorithm::M_Shake::m_shake_iteration
(
        topology::Topology const &topo,
        configuration::Configuration & conf,
        bool & convergence,
        int first,
        double dt,
        std::vector<topology::two_body_term_struct> const & constr,
        math::GenericMatrix<double> const & factor,
        math::GenericVec<double> const & constr_length2,
        math::GenericVec<math::Vec> const & dist_old
        ) {
  convergence = true;
  double const dt2 = dt * dt;

  // calculate constraints: solve A*l = c by direct matrix inversion

  // number of constraints -> only implemented for 3 constraints!!!
  unsigned int num_constr = constr.size();
  assert(num_constr == 3);

  // to store the positions and distances
  math::VArray &new_unconstr_pos = conf.current().pos;
  math::GenericVec<math::Vec> dist_new;
  math::GenericVec<double> dist2;
  math::GenericVec<double> diff;
  // to store matrix A
  math::GenericMatrix<double> A;

  // get the new distances and calculate the difference
  for (unsigned int k = 0; k < 3; ++k) {
    dist_new(k) = new_unconstr_pos(first+constr[k].i) - new_unconstr_pos(first+constr[k].j);
    dist2(k) = abs2(dist_new(k));
    diff(k) = constr_length2(k) - dist2(k);
  }

  // do we have to shake?
  if (fabs(diff(0)) >= constr_length2(0)*tolerance()*2.0 ||
      fabs(diff(1)) >= constr_length2(1)*tolerance()*2.0 ||
      fabs(diff(2)) >= constr_length2(2)*tolerance()*2.0) {

    // matrix A
    A(0,0) = dot(dist_old(0), dist_new(0))*factor(0,0);
    A(0,1) = dot(dist_old(1), dist_new(0))*factor(0,1);
    A(0,2) = dot(dist_old(2), dist_new(0))*factor(0,2);
    A(1,0) = dot(dist_old(0), dist_new(1))*factor(1,0);
    A(1,1) = dot(dist_old(1), dist_new(1))*factor(1,1);
    A(1,2) = dot(dist_old(2), dist_new(1))*factor(1,2);
    A(2,0) = dot(dist_old(0), dist_new(2))*factor(2,0);
    A(2,1) = dot(dist_old(1), dist_new(2))*factor(2,1);
    A(2,2) = dot(dist_old(2), dist_new(2))*factor(2,2);

    // vectors orthogonal? -> SHAKE error
    if (A(0,0) < constr_length2(0) * math::epsilon ||
        A(1,1) < constr_length2(1) * math::epsilon ||
        A(2,2) < constr_length2(2) * math::epsilon) {
      std::cout << "M_SHAKE ERROR\n"
                << "\tatom i    : " << first + 1 << "\n"
                << "\tatom j    : " << first + 2 << "\n"
                << "\tatom k    : " << first + 3 << "\n"
                << "\tfree i    : " << math::v2s(new_unconstr_pos(first)) << "\n"
                << "\tfree j    : " << math::v2s(new_unconstr_pos(first+1)) << "\n"
                << "\tfree k    : " << math::v2s(new_unconstr_pos(first+2)) << "\n"
                << "\tref r_ij     : " << math::v2s(dist_old(0)) << "\n"
                << "\tref r_ik     : " << math::v2s(dist_old(1)) << "\n"
                << "\tref r_jk     : " << math::v2s(dist_old(2)) << "\n"
                << "\tr_ij         : " << math::v2s(dist_new(0)) << "\n"
                << "\tr_ik         : " << math::v2s(dist_new(1)) << "\n"
                << "\tr_jk         : " << math::v2s(dist_new(2)) << "\n"
                << "\tA(i,i)        : " << A(0,0) << "\n"
                << "\tA(j,j)        : " << A(1,1) << "\n"
                << "\tA(k,k)        : " << A(2,2) << "\n"
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

    double f0 = (A(0,0)*diff(0) + A(0,1)*diff(1) + A(0,2)*diff(2)) / 2.0;
    double f1 = (A(1,0)*diff(0) + A(1,1)*diff(1) + A(1,2)*diff(2)) / 2.0;
    double f2 = (A(2,0)*diff(0) + A(2,1)*diff(1) + A(2,2)*diff(2)) / 2.0;

    DEBUG(10, "Lagrange multiplier of contraint 0: " << f0);
    DEBUG(10, "Lagrange multiplier of contraint 1: " << f1);
    DEBUG(10, "Lagrange multiplier of contraint 2: " << f2);

    math::Vec f01 = f0*dist_old(0) + f1*dist_old(1);
    math::Vec f02 = f2*dist_old(2) - f0*dist_old(0);
    math::Vec f12 = f1*dist_old(1) + f2*dist_old(2);

    // add to the forces
    conf.old().constraint_force(first) += f01;
    conf.old().constraint_force(first+1) += f02;
    conf.old().constraint_force(first+2) -= f12;

    // virial
    if (V == math::atomic_virial) {
        for (int a = 0; a < 3; ++a) {
          for (int aa = 0; aa < 3; ++aa) {
            conf.old().virial_tensor(a, aa) -=
                    (dist_old(0)(a) * dist_old(0)(aa) * f0 +
                     dist_old(1)(a) * dist_old(1)(aa) * f1 +
                     dist_old(2)(a) * dist_old(2)(aa) * f2) / dt2;
          }
        }
        DEBUG(12, "\tatomic virial done");
      }

    // calculate the new constrained positions
    new_unconstr_pos(first) += f01 / topo.mass()(first);
    new_unconstr_pos(first+1) += f02 / topo.mass()(first+1);
    new_unconstr_pos(first+2) -= f12 / topo.mass()(first+2);

    convergence = false;
  } // if we have to shake

  return 0;

}

/**
 * shake solvent.
 */
template<math::boundary_enum B, math::virial_enum V>
void algorithm::M_Shake
::solvent(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        double dt, int const max_iterations,
        int & error) {

  DEBUG(8, "\tshaking SOLVENT - M_SHAKE");

  if (!sim.mpi || m_rank == 0)
    m_timer.start("solvent");

  // the first atom of a solvent
  unsigned int first = topo.num_solute_atoms();

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int tot_iterations = 0;

  error = 0;
  int my_error = error;

  const unsigned int num_atoms = topo.num_atoms();
  math::Periodicity<B> periodicity(conf.current().box);

  /*#ifdef XXMPI
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
  #endif*/

  // for all solvents
  for (unsigned int i = 0; i < topo.num_solvents(); ++i) {
    const unsigned int num_solvent_atoms = topo.solvent(i).num_atoms();

    math::GenericMatrix<double> factor;
    unsigned int k = 0;
    unsigned int num_constr = topo.solvent(i).distance_constraints().size();
    math::GenericVec<double> constr_length2;
    math::GenericVec<math::Vec> dist_old;

    
    for (typename std::vector<topology::two_body_term_struct>::const_iterator
      it_k = topo.solvent(i).distance_constraints().begin(),
      to_k = topo.solvent(i).distance_constraints().end(); it_k != to_k; ++it_k, ++k) {

      assert(parameter().size() > it_k->type);
      constr_length2(k) = parameter()[it_k->type].r0 * parameter()[it_k->type].r0;
      unsigned int l = 0;
      for (typename std::vector<topology::two_body_term_struct>::const_iterator
        it_l = topo.solvent(i).distance_constraints().begin(),
        to_l = topo.solvent(i).distance_constraints().end(); it_l != to_l; ++it_l, ++l) {
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
        factor(k, l) = ((d11 - d12) / topo.mass()(first+it_k->i) + (d22 - d21) /
                                      topo.mass()(first+it_k->j));
      }
    }

    // loop over the molecules
    for (unsigned int nm = 0; nm < topo.num_solvent_molecules(i);
            ++nm, first += num_solvent_atoms) {

      /*#ifdef XXMPI
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
      #endif*/

      k = 0;
      for (typename std::vector<topology::two_body_term_struct>::const_iterator
      it_k = topo.solvent(i).distance_constraints().begin(),
      to_k = topo.solvent(i).distance_constraints().end(); it_k != to_k; ++it_k, ++k) {
        dist_old(k) = conf.old().pos(first+it_k->i) - conf.old().pos(first+it_k->j);
      }

      int num_iterations = 0;
      bool convergence = false;
      while (!convergence) {
        DEBUG(9, "\titeration" << std::setw(10) << num_iterations
                << " ma_iterations" << std::setw(10) << max_iterations);

        if (m_shake_iteration<B, V >
                (topo, conf, convergence, first, dt,
                topo.solvent(i).distance_constraints(),
                factor, constr_length2, dist_old)) {

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

    } // molecules
    if (my_error != error) break;

  } // solvents

  // reduce everything
  /*#ifdef XXMPI
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
  #endif*/
  error = my_error;

  // constraint force
  const double dt2 = sim.time_step_size() * sim.time_step_size();
  for (unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
    conf.old().constraint_force(i) *= 1 / dt2;
    DEBUG(5, "constraint_force " << math::v2s(conf.old().constraint_force(i)));
  }

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
  std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
          to = constrained_atoms().end();
  for (; it != to; ++it) {
    conf.old().constraint_force(*it) = 0.0;
  }

  math::SArray masses;

  /*if (sim.param().posrest.posrest == simulation::posrest_const) {
    masses = topo.mass();
    std::vector<topology::position_restraint_struct>::const_iterator
    it = topo.position_restraints().begin(),
            to = topo.position_restraints().end();
    for (; it != to; ++it) {
      topo.mass()[it->seq] = std::numeric_limits<double>::infinity();
    }
  }*/

  // check whether we shake

  if (sim.param().system.nsm &&
          sim.param().constraint.solvent.algorithm == simulation::constr_m_shake) {

    DEBUG(8, "\twe need to shake SOLVENT");
    SPLIT_VIRIAL_BOUNDARY(solvent,
            topo, conf, sim, sim.time_step_size(),
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
    std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.current().vel(*it) = (conf.current().pos(*it) - conf.old().pos(*it)) /
              sim.time_step_size();
    }
  }

  /*if (sim.param().posrest.posrest == simulation::posrest_const) {
    std::vector<topology::position_restraint_struct>::const_iterator
    it = topo.position_restraints().begin(),
            to = topo.position_restraints().end();
    for (; it != to; ++it) {
      topo.mass()[it->seq] = masses[it->seq];
    }
  }*/

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

  /*#ifdef XXMPI
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
  #endif*/

  if (sim.param().constraint.solute.algorithm == simulation::constr_m_shake) {
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
  }

  if (sim.param().constraint.solvent.algorithm == simulation::constr_m_shake) {
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
            "m_shake", io::message::error);
  }

  if (!quiet)
    os << "END\n";

  return 0;
}
