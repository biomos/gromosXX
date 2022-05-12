/**
 * @file gpu_shake.cc
 * contains the methods for
 * the class GPU_Shake.
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

#include "../../algorithm/constraints/gpu_shake.h"
#include "../../algorithm/constraints/gpu_shake_thread.h"
#ifdef HAVE_LIBCUDART
#include <cudaKernel.h>
#endif

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"
#include <cassert>

#include <limits>
#include <vector>

#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

#define GPU_ID 0
/**
 * Constructor.
 */
algorithm::GPU_Shake
::GPU_Shake(double const tolerance, int const max_iterations,
        std::string const name)
: Algorithm(name),
m_tolerance(tolerance),
m_max_iterations(max_iterations) {
}

/**
 * Destructor.
 */
algorithm::GPU_Shake
::~GPU_Shake() {
}

void algorithm::GPU_Shake
::tolerance(double const tol) {
  m_tolerance = tol;
}

/**
 * shake solvent.
 */
void algorithm::GPU_Shake
::solvent(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        double dt, int const max_iterations,
        int & error) {

  DEBUG(8, "\tshaking SOLVENT - GPU_SHAKE");

  if (!sim.mpi || m_rank == 0)
    m_timer.start("solvent");

  error = 0;

#ifdef OMP
  const unsigned int num_gpus = sim.param().constraint.solvent.number_gpus;
  DEBUG (10, "GPU_SHAKE : Enter OMP")
#pragma omp parallel
  {
    unsigned int tid = omp_get_thread_num();
    if (tid < num_gpus) {   // then shake!
      int shake_fail_mol = -1;
      DEBUG(10, "GPU_SHAKE : OMP " << tid << " : Apply")
      shake_fail_mol = m_shake_set[tid]->apply();

      // if error
      DEBUG(10, "GPU_Shake : Is there a fail molecule?")
      if (shake_fail_mol > 0) {
        const unsigned int first = topo.num_solute_atoms() + (shake_fail_mol - 1) * topo.solvent(0).num_atoms();
        std::cout << "GPU_SHAKE ERROR\n"
                << "\tatom i    : " << first + 1 << "\n"
                << "\tatom j    : " << first + 2 << "\n"
                << "\tatom k    : " << first + 3 << "\n";
        error = 1;
      } // error handling

    } // for all gpus
  } // OMP
  DEBUG(10, "GPU_SHAKE : OMP Finished");

#else
  std::cerr << "OMP not defined, why are we here???" << std::endl;
#endif

  if (!sim.mpi || m_rank == 0)
    m_timer.stop("solvent");
} // shake solvent


/**
 * apply the GPU_SHAKE algorithm
 */
int algorithm::GPU_Shake::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(7, "applying GPU_SHAKE");
  if (!sim.mpi || m_rank == 0)
    m_timer.start();

  int error = 0;

  // check whether we shake

  if (sim.param().system.nsm &&
          sim.param().constraint.solvent.algorithm == simulation::constr_gpu_shake) {

    DEBUG(8, "\twe need to shake SOLVENT");
    solvent(topo, conf, sim, sim.time_step_size(),
            m_max_iterations, error);
    if (error) {
      std::cout << "GPU_SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLVENT "
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
    const unsigned int num_atoms = topo.num_atoms();
    const double dti = 1.0 / sim.time_step_size();
    for (unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
      conf.current().vel(i) = (conf.current().pos(i) - conf.old().pos(i)) * dti;
    }
  }

  if (!sim.mpi || m_rank == 0)
    m_timer.stop();
  // return success!
  return 0;

}

int algorithm::GPU_Shake::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) {
  
  // Make constraint matrix
  // the first atom of a solvent
  const unsigned int first = topo.num_solute_atoms();

  math::GenericMatrix<double> factor;
  math::GenericVec<double> constr_length2;

  unsigned int k = 0;
  for (std::vector<topology::two_body_term_struct>::const_iterator
    it_k = topo.solvent(0).distance_constraints().begin(),
          to_k = topo.solvent(0).distance_constraints().end(); it_k != to_k; ++it_k, ++k) {

    assert(topo.bond_types_harm().size() > it_k->type);
    constr_length2(k) = topo.bond_types_harm()[it_k->type].r0 * topo.bond_types_harm()[it_k->type].r0;
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
      factor(k, l) = ((d11 - d12) / topo.mass()(first + it_k->i) + (d22 - d21) /
              topo.mass()(first + it_k->j));
    }
  }

  // Make the mass
  math::Vec mass(topo.mass()(first), topo.mass()(first + 1), topo.mass()(first + 2));

#ifdef OMP
  const unsigned int num_gpus = sim.param().constraint.solvent.number_gpus;
#pragma omp parallel
  {
    unsigned int tid = omp_get_thread_num();
    if (tid < num_gpus) {     // then create a thread for the gpu
      GPU_Shake_Thread * gpu_shake = nullptr;
      gpu_shake = new GPU_Shake_Thread(topo, conf, sim, factor, mass, constr_length2, tolerance(), num_gpus, tid, os, quiet);
      if (gpu_shake == NULL)
        std::cerr << "Could not allocate space for GPU_Shake_Thread!" << std::endl;

#pragma omp critical
      {
        DEBUG(5, "GPU_Shake : Pushed back GPU_Thread " << tid)
        m_shake_set.push_back(gpu_shake);
      }
    } // for all gpus

  } // OMP
#else
  std::cerr << "OMP not defined, why are we here???" << std::endl;
  return E_ILLEGAL;
#endif  

    if (!quiet) {
    os << "GPU-SHAKE\n"
            << "\tsolvent\t";

    if (sim.param().constraint.solvent.algorithm == simulation::constr_gpu_shake) {
      if (sim.mpi)
        os << "ON (MPI parallel version)\n";
      else
        os << "ON\n";
      os << "\t\ttolerance = "
              << sim.param().constraint.solvent.shake_tolerance << "\n";
    } else os << "OFF\n";
  }


  m_rank = 0;
  m_size = 1;

  if (sim.param().start.shake_pos) {
    if (!quiet)
      os << "\n\tshaking initial positions\n";

    // old and current pos and vel are the same for constrained atoms...

    const unsigned int num_atoms = topo.num_atoms();
    for (unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
      conf.old().pos(i) = conf.current().pos(i);
      conf.old().vel(i) = conf.current().vel(i);
    }

    // shake the current ones
    if (apply(topo, conf, sim))
      return E_SHAKE_FAILURE;

    for (unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
      // restore the velocities
      conf.current().vel(i) = conf.old().vel(i);
      // take a step back
      conf.old().pos(i) = conf.current().pos(i);
    }

    if (sim.param().start.shake_vel) {
      if (!quiet)
        os << "\tshaking initial velocities\n";

      const unsigned int num_atoms = topo.num_atoms();
      for (unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
        conf.current().pos(i) = conf.old().pos(i) -
                sim.time_step_size() * conf.old().vel(i);
      }

      // shake again
      if (apply(topo, conf, sim))
        return E_SHAKE_FAILURE;

      for (unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
        // restore the positions
        conf.current().pos(i) = conf.old().pos(i);
        // velocities are in opposite direction (in time)
        conf.current().vel(i) = -1.0 * conf.current().vel(i);
        conf.old().vel(i) = conf.current().vel(i);
      }
    } // if shake vel
  } else if (sim.param().start.shake_vel) {
    io::messages.add("shaking velocities without shaking positions illegal.",
            "GPU_Shake", io::message::error);
  }

  // set the constraint force to zero because it is not calculated and uninitialized values are confusing
  const unsigned int num_atoms = topo.num_atoms();
  for (unsigned int i = topo.num_solute_atoms(); i < num_atoms; ++i) {
    conf.current().constraint_force(i) = 0.0;
    conf.old().constraint_force(i) = 0.0;
  }

  if (!quiet)
    os << "END\n";

  return 0;
}

