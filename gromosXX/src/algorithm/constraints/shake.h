/**
 * @file shake.h
 * the shake algorithm.
 */

#ifndef INCLUDED_SHAKE_H
#define INCLUDED_SHAKE_H

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm
{
  /**
   * @class Shake
   * implements the shake algorithm.
   */
  class Shake : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Shake(double const tolerance = 0.000001, 
	  int const max_iterations = 1000,
	  std::string const name = "Shake");

    /**
     * Destructor.
     */
    virtual ~Shake();
        
    /**
     * apply shake.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

    /**
     * tolerance.
     */
    double const & tolerance()const {return m_tolerance;}
    /**
     * max iterations.
     */
    int const & max_iterations()const {return m_max_iterations;}

    /**
     * the const bond type parameter.
     */
    std::vector<interaction::bond_type_struct> const &parameter()const
    {
      return m_parameter;
    }
    /**
     * the bond type parameter.
     */
    std::vector<interaction::bond_type_struct> & parameter()
    {
      return m_parameter;
    }
    /**
     * accessor to the constrained atoms
     */
    std::set<unsigned int> & constrained_atoms() {
      return m_constrained_atoms;
    }
     /**
     * accessor to the constrained atoms
     */
    const std::set<unsigned int> & constrained_atoms() const {
      return m_constrained_atoms;
    }

    /**
     * initialize startup positions and velocities
     * if required.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

  protected:

    /**
     * shake tolerance
     */
    double m_tolerance;
    /**
     * max iterations
     */
    const int m_max_iterations;
    /**
     * bond parameter
     */
    std::vector<interaction::bond_type_struct> m_parameter;
    /**
     * the atoms that are involved in the contraints
     */
    std::set<unsigned int> m_constrained_atoms;
    /** 
     * rank and size for parallelization
     */
    int m_rank, m_size;

    template<math::boundary_enum B, math::virial_enum V>
    int shake_iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     bool & convergence,
     int first,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     std::vector<topology::two_body_term_struct> const & constr,
     double dt,
     math::Periodicity<B> const & periodicity
     );
    
    
    template<math::boundary_enum B, math::virial_enum V>
    int dih_constr_iteration
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim,
     bool & convergence,
     std::vector<bool> & skip_now,
     std::vector<bool> & skip_next,
     math::Periodicity<B> const & periodicity
     );


    template<math::boundary_enum B, math::virial_enum V>
    void solute
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim,
     int const max_iterations,
     int & error
     );
    
    template<math::boundary_enum B, math::virial_enum V>
    void solvent
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     double dt, int const max_iterations,
     int & error
     );

  };
  
} //algorithm



/*********************************************************************
 * Definitions of templated functions
 *********************************************************************/

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

// save previous values
#ifdef MODULE
#define MODULE_PREV MODULE
#endif

#ifdef SUBMODULE
#define SUBMODULE_PREV SUBMODULE
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

// dihedral constraints template method
#include "dihedral_constraint.cc"


/**
 * do one iteration
 */
template<math::boundary_enum B, math::virial_enum V>
int algorithm::Shake::shake_iteration
(
        topology::Topology const &topo,
        configuration::Configuration & conf,
        bool & convergence,
        int first,
        std::vector<bool> &skip_now,
        std::vector<bool> &skip_next,
        std::vector<topology::two_body_term_struct> const & constr,
        double dt,
        math::Periodicity<B> const & periodicity
        )
 {
  convergence = true;

  // index for constraint_force...
  unsigned int k = 0;
  double const dt2 = dt * dt;

  // and constraints
  for (typename std::vector<topology::two_body_term_struct>
          ::const_iterator
          it = constr.begin(),
          to = constr.end();
          it != to;
          ++it, ++k) {

    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j]) continue;
    if (topo.inverse_mass(first + it->i) == 0 &&
            topo.inverse_mass(first + it->j) == 0) continue;

    DEBUG(10, "i: " << it->i << " j: " << it->j << " first: " << first);

    // the position
    math::Vec &pos_i = conf.current().pos(first + it->i);
    math::Vec &pos_j = conf.current().pos(first + it->j);

    DEBUG(10, "\ni: " << math::v2s(pos_i) << "\nj: " << math::v2s(pos_j));

    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);
    DEBUG(12, "ni:  " << math::v2s(r));

    double dist2 = abs2(r);

    assert(parameter().size() > it->type);
    double constr_length2 = parameter()[it->type].r0 * parameter()[it->type].r0;
    double diff = constr_length2 - dist2;

    DEBUG(13, "constr: " << constr_length2 << " dist2: " << dist2);

    if (fabs(diff) >= constr_length2 * tolerance() * 2.0) {
      // we have to shake
      DEBUG(10, "shaking");

      // the reference position
      const unsigned int atom_i = first + it->i;
      const unsigned int atom_j = first + it->j;
      const math::Vec &ref_i = conf.old().pos(atom_i);
      const math::Vec &ref_j = conf.old().pos(atom_j);

      math::Vec ref_r;
      periodicity.nearest_image(ref_i, ref_j, ref_r);

      double sp = dot(ref_r, r);

      DEBUG(5, "ref i " << math::v2s(ref_i) << " ref j " << math::v2s(ref_j));
      DEBUG(5, "free i " << math::v2s(pos_i) << " free j " << math::v2s(pos_j));
      DEBUG(5, "ref r " << math::v2s(ref_r));
      DEBUG(5, "r " << math::v2s(r));

      if (sp < constr_length2 * math::epsilon) {
        /*
        io::messages.add("SHAKE error. vectors orthogonal",
             "Shake::???",
             io::message::critical);
         */
        std::cout << "SHAKE ERROR\n"
                << "\tatom i    : " << atom_i + 1 << "\n"
                << "\tatom j    : " << atom_j + 1 << "\n"
                // << "\tfirst     : " << first << "\n"
                << "\tref i     : " << math::v2s(ref_i) << "\n"
                << "\tref j     : " << math::v2s(ref_j) << "\n"
                << "\tfree i    : " << math::v2s(pos_i) << "\n"
                << "\tfree j    : " << math::v2s(pos_j) << "\n"
                << "\tref r     : " << math::v2s(ref_r) << "\n"
                << "\tr         : " << math::v2s(r) << "\n"
                << "\tsp        : " << sp << "\n"
                << "\tconstr    : " << constr_length2 << "\n"
                << "\tdiff      : " << diff << "\n"
                << "\tforce i   : " << math::v2s(conf.old().force(atom_i)) << "\n"
                << "\tforce j   : " << math::v2s(conf.old().force(atom_j)) << "\n"
                << "\tvel i     : " << math::v2s(conf.current().vel(atom_i)) << "\n"
                << "\tvel j     : " << math::v2s(conf.current().vel(atom_j)) << "\n"
                << "\told vel i : " << math::v2s(conf.old().vel(atom_i)) << "\n"
                << "\told vel j : " << math::v2s(conf.old().vel(atom_j)) << "\n\n";

        return E_SHAKE_FAILURE;
      }

      // lagrange multiplier
      double lambda = diff / (sp * 2.0 *
              (topo.inverse_mass()(atom_i) +
               topo.inverse_mass()(atom_j)));

      DEBUG(10, "lagrange multiplier " << lambda);

      const math::Vec cons_force = lambda * ref_r;
      conf.old().constraint_force(atom_i) += cons_force;
      conf.old().constraint_force(atom_j) -= cons_force;

      if (V == math::atomic_virial) {
        for (int a = 0; a < 3; ++a) {
          for (int aa = 0; aa < 3; ++aa) {
            conf.old().virial_tensor(a, aa) +=
                    ref_r(a) * ref_r(aa) * lambda / dt2;
          }
        }
        DEBUG(12, "\tatomic virial done");
      }

      // update positions
      ref_r *= lambda;
      pos_i += ref_r * topo.inverse_mass()(first + it->i);
      pos_j -= ref_r * topo.inverse_mass()(first + it->j);

      convergence = false;

      // consider atoms in the next step
      skip_next[it->i] = false;
      skip_next[it->j] = false;

    } // we have to shake
  } // constraints


  return 0;

}

/**
 * shake solute
 */
template<math::boundary_enum B, math::virial_enum V>
void algorithm::Shake::
solute(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation const & sim,
        int const max_iterations,
        int & error) {
  // for now shake the whole solute in one go,
  // not bothering about submolecules...

  DEBUG(8, "\tshaking SOLUTE");
  math::Periodicity<B> periodicity(conf.current().box);

  m_timer.start("solute");

  const unsigned int num_atoms = topo.num_solute_atoms();
  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  int first = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);

  bool convergence = false;
  while (!convergence) {
    DEBUG(9, "\titeration" << std::setw(10) << num_iterations);

    // distance constraints
    bool dist_convergence = true;

    if (topo.solute().distance_constraints().size() &&
            sim.param().constraint.solute.algorithm == simulation::constr_shake &&
            sim.param().constraint.ntc > 1) {

      DEBUG(7, "SHAKE: distance constraints iteration");

      if (shake_iteration<B, V >
              (topo, conf, dist_convergence, first, skip_now, skip_next,
              topo.solute().distance_constraints(), sim.time_step_size(),
              periodicity)
              ) {
        io::messages.add("SHAKE error. vectors orthogonal",
                "Shake::solute",
                io::message::error);
        std::cout << "SHAKE failure in solute!" << std::endl;
        error = E_SHAKE_FAILURE_SOLUTE;
        return;
      }
    }

    // dihedral constraints
    bool dih_convergence = true;
    if (sim.param().dihrest.dihrest == simulation::dihedral_constr) {

      DEBUG(7, "SHAKE: dihedral constraints iteration");

      if (dih_constr_iteration<B, V >
              (topo, conf, sim, dih_convergence, skip_now, skip_next, periodicity)
              ) {
        io::messages.add("SHAKE error: dihedral constraints",
                "Shake::solute",
                io::message::error);
        std::cout << "SHAKE failure in solute dihedral constraints!" << std::endl;
        error = E_SHAKE_FAILURE_SOLUTE;
        return;
      }
    }

    convergence = dist_convergence && dih_convergence;

    if (++num_iterations > max_iterations) {
      io::messages.add("SHAKE error. too many iterations",
              "Shake::solute",
              io::message::error);
      error = E_SHAKE_FAILURE_SOLUTE;
      return;
    }

    skip_now = skip_next;
    skip_next.assign(skip_next.size(), true);

  } // convergence?

  // constraint force
  const double dt2 = sim.time_step_size() * sim.time_step_size();
  for (unsigned int i = 0; i < num_atoms; ++i) {
    conf.old().constraint_force(i) *= 1 / dt2;
    DEBUG(5, "constraint_force " << math::v2s(conf.old().constraint_force(i)));
  }
  error = 0;

  m_timer.stop("solute");

} // solute
/**
 * shake solvent.
 */
template<math::boundary_enum B, math::virial_enum V>
void algorithm::Shake
::solvent(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        double dt, int const max_iterations,
        int & error) {

  DEBUG(8, "\tshaking SOLVENT");

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
  for (unsigned int i = 0; i < topo.num_solvents(); ++i) {
    const unsigned int num_solvent_atoms = topo.solvent(i).num_atoms();
    // loop over the molecules
    for (unsigned int nm = 0; nm < topo.num_solvent_molecules(i);
            ++nm, first += num_solvent_atoms) {

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

      skip_now.assign(num_solvent_atoms, false);
      skip_next.assign(num_solvent_atoms, true);

      int num_iterations = 0;
      bool convergence = false;
      while (!convergence) {
        DEBUG(9, "\titeration" << std::setw(10) << num_iterations);

        if (shake_iteration<B, V >
                (topo, conf, convergence, first, skip_now, skip_next,
                topo.solvent(i).distance_constraints(), dt,
                periodicity)) {

          io::messages.add("SHAKE error. vectors orthogonal",
                  "Shake::solvent", io::message::error);

          std::cout << "SHAKE failure in solvent!" << std::endl;
          my_error = E_SHAKE_FAILURE_SOLVENT;
          break;
        }

        // std::cout << num_iterations+1 << std::endl;
        if (++num_iterations > max_iterations) {
          io::messages.add("SHAKE error. too many iterations",
                  "Shake::solvent",
                  io::message::critical);
          my_error = E_SHAKE_FAILURE_SOLVENT;
          break;
        }

        skip_now = skip_next;
        skip_next.assign(skip_next.size(), true);

      } // while(!convergence)
      if (my_error != error) break;

      tot_iterations += num_iterations;

    } // molecules
    if (my_error != error) break;

  } // solvents

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


// set macros to previous values if needed
#ifdef MODULE_PREV
#undef MODULE
#define MODULE MODULE_PREV
#undef MODULE_PREV
#endif

#ifdef SUBMODULE_PREV
#undef SUBMODULE
#define SUBMODULE SUBMODULE_PREV
#undef SUBMODULE_PREV
#endif

#endif
