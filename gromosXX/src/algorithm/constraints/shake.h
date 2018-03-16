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
   * @struct ConstraintGroup
   * grouped constraints
   */
  struct ConstraintGroup {
    std::vector<topology::dihedral_restraint_struct> dihedral_restraints;
    std::vector<topology::two_body_term_struct> distance_restraints;
#ifdef XXMPI
    std::vector<unsigned int> unaffected_indices;
#endif
  };

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
     * the atoms that are involved in the contraints
     */
    std::set<unsigned int> m_constrained_atoms;

    /**
     * disjoined constrained groups
     */
    std::vector<ConstraintGroup> m_constraint_groups;
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
     std::vector<topology::dihedral_restraint_struct> const & dihedral_restraints,
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

    assert(topo.bond_types_harm().size() > it->type);
    double constr_length2 = topo.bond_types_harm()[it->type].r0 * topo.bond_types_harm()[it->type].r0;
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

      DEBUG(10, "ref i " << math::v2s(ref_i) << " ref j " << math::v2s(ref_j));
      DEBUG(10, "free i " << math::v2s(pos_i) << " free j " << math::v2s(pos_j));
      DEBUG(10, "ref r " << math::v2s(ref_r));
      DEBUG(10, "r " << math::v2s(r));

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

  if (!sim.mpi || m_rank == 0)
    m_timer.start("solute");

  const unsigned int num_atoms = topo.num_solute_atoms();
  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  skip_next.assign(topo.solute().num_atoms(), true);
  skip_now.assign(topo.solute().num_atoms(), false);

  int first = 0;
  error = 0;
  int my_error = error;

  conf.special().constraints.num_solute_dist_iterations=0;
  conf.special().constraints.num_dih_iterations=0;

  const unsigned int group_id = m_rank;
  int num_iterations = 0;
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
              m_constraint_groups[group_id].distance_restraints, sim.time_step_size(),
              periodicity)
              ) {
        io::messages.add("SHAKE error. vectors orthogonal",
                "Shake::solute",
                io::message::error);
        std::cout << "SHAKE failure in solute!" << std::endl;
        my_error = E_SHAKE_FAILURE_SOLUTE;
        break;
      }
      if (!dist_convergence) conf.special().constraints.num_solute_dist_iterations++;
    }

    // dihedral constraints
    bool dih_convergence = true;
    if (sim.param().dihrest.dihrest == simulation::dihedral_constr) {

      DEBUG(7, "SHAKE: dihedral constraints iteration");

      if (dih_constr_iteration<B, V >
              (topo, conf, sim, dih_convergence, skip_now, skip_next, m_constraint_groups[group_id].dihedral_restraints, periodicity)
              ) {
        io::messages.add("SHAKE error: dihedral constraints",
                "Shake::solute",
                io::message::error);
        std::cout << "SHAKE failure in solute dihedral constraints!" << std::endl;
        my_error = E_SHAKE_FAILURE_SOLUTE;
        break;
      }
    }
    if (!dih_convergence) conf.special().constraints.num_dih_iterations++;

    convergence = dist_convergence && dih_convergence;

    if (++num_iterations > max_iterations) {
      io::messages.add("SHAKE error. too many iterations",
              "Shake::solute",
              io::message::error);
      my_error = E_SHAKE_FAILURE_SOLUTE;
      break;
    }
    conf.special().constraints.num_iterations = num_iterations;
    std::swap(skip_next, skip_now);
    skip_next.assign(skip_next.size(), true);
  } // convergence?
  DEBUG(10, "SHAKE: num_iterations " << num_iterations << " num_dih_iterations " << conf.special().constraints.num_dih_iterations);

  // reduce errors
#ifdef XXMPI
  if (sim.mpi) {
    MPI::COMM_WORLD.Allreduce(&my_error, &error, 1, MPI::INT, MPI::MAX);
  } else
  error = my_error;
#else
  error = my_error;
#endif

  //TODO: this is just temporary to write out constrained dih.angles to trs
  // if we want to keep it we should do it in a more efficient way probably
  if (conf.special().dihedralres.d.size() > 0) {
  math::VArray &pos   = conf.current().pos;

  std::vector<topology::dihedral_restraint_struct>::const_iterator
    it = m_constraint_groups[group_id].dihedral_restraints.begin(),
    to = m_constraint_groups[group_id].dihedral_restraints.end();

  for( int i=0; it != to; ++it, ++i){
    // calculate dihedral angle
    DEBUG(9, "dihedral angle " << it->i << "-" << it->j << "-" << it->k << "-" << it->l);

    math::Vec r12, r32, r34;
    periodicity.nearest_image(pos(it->i), pos(it->j), r12);
    periodicity.nearest_image(pos(it->k), pos(it->j), r32);
    periodicity.nearest_image(pos(it->k), pos(it->l), r34);

    // eq 35
    const math::Vec r52 = math::cross(r12, r32);
    const double d52 = math::abs(r52);
    // eq 36
    const math::Vec r63 = math::cross(r32, r34);
    const double d63 = math::abs(r63);
    // eq 37
    const int sign_phi = (math::dot(r12, r63) >= 0.0) ? 1 : -1;
    // eq 34
    const double cos_phi = math::dot(r52, r63) / (d52 * d63);

    double phi;
    // cos_phi can be >1 or <-1 because of precision limits
    if (cos_phi > 1) phi=0.0;
    else if (cos_phi < -1) phi=math::Pi;
    else phi = sign_phi * acos(cos_phi);

    while(phi < it->phi - math::Pi)
      phi += 2 * math::Pi;
    while(phi > it->phi + math::Pi)
      phi -= 2 * math::Pi;
      
    conf.special().dihedralres.d[i] = phi;
    conf.special().dihedralres.energy[i] = 0.0;
  }
  }

  // constraint force
  const double dt2 = sim.time_step_size() * sim.time_step_size();
  for (unsigned int i = 0; i < num_atoms; ++i) {
    conf.old().constraint_force(i) *= 1 / dt2;
    DEBUG(5, "constraint_force " << math::v2s(conf.old().constraint_force(i)));
  }
  // error = 0;

  if (!sim.mpi || m_rank == 0)
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

  conf.special().constraints.num_solvent_dist_iterations = 0;

  const unsigned int num_atoms = topo.num_atoms();
  math::Periodicity<B> periodicity(conf.current().box);

  // for all solvents
  for (unsigned int i = 0; i < topo.num_solvents(); ++i) {
    const unsigned int num_solvent_atoms = topo.solvent(i).num_atoms();
    // loop over the molecules
    for (unsigned int nm = 0; nm < topo.num_solvent_molecules(i);
            ++nm, first += num_solvent_atoms) {

#ifdef XXMPI
      math::VArray & pos = conf.current().pos;
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
      conf.special().constraints.num_solvent_dist_iterations=tot_iterations;

    } // molecules
    if (my_error != error) break;

  } // solvents

  // reduce errors
#ifdef XXMPI
  if (sim.mpi) {
    if (m_rank == 0) {
      // Master 
      // reduce the error to all processors
      MPI::COMM_WORLD.Allreduce(&my_error, &error, 1, MPI::INT, MPI::MAX);
    } else {
      // reduce the error to all processors
      MPI::COMM_WORLD.Allreduce(&my_error, &error, 1, MPI::INT, MPI::MAX);
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
