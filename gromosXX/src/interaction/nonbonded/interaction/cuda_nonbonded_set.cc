/**
 * @file cuda_nonbonded_set.cc
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>
#include <interaction/nonbonded/interaction/cuda_nonbonded_set.h>
#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>

#include <interaction/nonbonded/interaction/latticesum.h>

#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>

#include <math/periodicity.h>
#include <math/volume.h>

#include <interaction/nonbonded/interaction/nonbonded_set.h>

#include <util/debug.h>

#include <configuration/energy.h>

#include "storage.h"

#ifdef HAVE_LIBCUKERNEL
#include <cudaKernel.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::CUDA_Nonbonded_Set::~CUDA_Nonbonded_Set() {
#ifdef HAVE_LIBCUKERNEL
  cudakernel::CleanUp();
#endif
}
interaction::CUDA_Nonbonded_Set
::CUDA_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param,
		int rank, int num_threads)
  : Nonbonded_Set(pairlist_alg, param, rank, num_threads) ,
    m_outerloop(param) {

  // these values should be fine for a cutoff of 0.8 / 1.4. They have
  // to be increased for larger cutoffs. 
  // unfortunately there is no way to check these values at runtime so you
  // have to make sure youself that they are big enough.
  estNeigh_long = 500;
  estNeigh_short = 400;
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::CUDA_Nonbonded_Set
::calculate_interactions(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(4, "Nonbonded_Set::calculate_interactions");

  int error = 0;

  m_storage.zero();
  m_storage_cuda.zero();
  const bool pairlist_update = !(sim.steps() % sim.param().pairlist.skip_step);


  if (pairlist_update) {
    DEBUG(6, "\tdoing longrange...");

    m_longrange_storage.zero();
    m_longrange_storage_cuda.zero();
#ifdef HAVE_LIBCUKERNEL
    m_pairlist_alg.timer().start("pairlist cuda");
    double * Pos = &conf.current().pos(topo.num_solute_atoms())(0);
    cudakernel::cudaCalcPairlist(Pos);
    m_pairlist_alg.timer().stop("pairlist cuda");
#endif

    m_pairlist_alg.update(topo, conf, sim,
            pairlist(),
            m_rank, topo.num_atoms(), m_num_threads);

  }

  if (pairlist_update) {
#ifdef HAVE_LIBCUKERNEL
    const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
#endif
    DEBUG(8, "doing longrange calculation");

    // in case of LS only LJ are calculated
    m_pairlist_alg.timer().start("longrange");

    m_outerloop.lj_crf_outerloop(topo, conf, sim,
            m_pairlist.solute_long, m_pairlist.solvent_long,
            m_longrange_storage);

    m_pairlist_alg.timer().stop("longrange");

#ifdef HAVE_LIBCUKERNEL
    m_pairlist_alg.timer().start("longrange-cuda");
    double * For = &m_longrange_storage.force(topo.num_solute_atoms())(0);
    double * e_lj = &m_longrange_storage.energies.lj_energy[egroup][egroup];
    double * e_crf = &m_longrange_storage.energies.crf_energy[egroup][egroup];
    double * Pos = &conf.current().pos(topo.num_solute_atoms())(0);
    error += cudakernel::cudaCalcForces(Pos, For, e_lj, e_crf, true);
    m_pairlist_alg.timer().stop("longrange-cuda");
#endif
  }
  // calculate forces / energies
  DEBUG(6, "\tshort range interactions");

  m_pairlist_alg.timer().start("shortrange");

  m_outerloop.lj_crf_outerloop(topo, conf, sim,
          m_pairlist.solute_short, m_pairlist.solvent_short,
          m_storage);

  m_pairlist_alg.timer().stop("shortrange");

#ifdef HAVE_LIBCUKERNEL

  m_pairlist_alg.timer().start("shortrange-cuda");
  double * Pos = &conf.current().pos(topo.num_solute_atoms())(0);
  double * For = &m_storage.force(topo.num_solute_atoms())(0);
  const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
  double * e_lj = &m_storage.energies.lj_energy[egroup][egroup];
  double * e_crf = &m_storage.energies.crf_energy[egroup][egroup];
  cudakernel::cudaCalcForces(Pos, For, e_lj, e_crf, false);
  m_pairlist_alg.timer().stop("shortrange-cuda");
#endif

  if (m_rank == 0 && error > 0)
    return 1;

  DEBUG(6, "\t1,4 - interactions");
  m_pairlist_alg.timer().start("1,4 interaction");
  m_outerloop.one_four_outerloop(topo, conf, sim, m_storage);
  m_pairlist_alg.timer().stop("1,4 interaction");

  // possibly do the RF contributions due to excluded atoms
  if (sim.param().nonbonded.rf_excluded) {
    std::cout << "rf_excluded" << std::endl;
    DEBUG(7, "\tRF excluded interactions and self term");
    m_pairlist_alg.timer().start("RF excluded interaction");
    m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_storage);
    m_pairlist_alg.timer().stop("RF excluded interaction");
  }


  // add long-range force
  DEBUG(6, "\t(set) add long range forces");

  m_storage.force += m_longrange_storage.force;
  m_storage_cuda.force += m_longrange_storage_cuda.force;

  // and long-range energies
  DEBUG(6, "\t(set) add long range energies");
  const unsigned int lj_e_size = unsigned(m_storage.energies.lj_energy.size());

  for (unsigned int i = 0; i < lj_e_size; ++i) {
    for (unsigned int j = 0; j < lj_e_size; ++j) {
      m_storage.energies.lj_energy[i][j] +=
              m_longrange_storage.energies.lj_energy[i][j];
      m_storage.energies.crf_energy[i][j] +=
              m_longrange_storage.energies.crf_energy[i][j];
    }
  }

  for (unsigned int i = 0; i < lj_e_size; ++i) {
    for (unsigned int j = 0; j < lj_e_size; ++j) {
      m_storage_cuda.energies.lj_energy[i][j] +=
              m_longrange_storage_cuda.energies.lj_energy[i][j];
      m_storage_cuda.energies.crf_energy[i][j] +=
              m_longrange_storage_cuda.energies.crf_energy[i][j];
    }
  }

  // add longrange virial
  if (sim.param().pcouple.virial) {
    DEBUG(6, "\t(set) add long range virial");

    m_storage.virial_tensor += m_longrange_storage.virial_tensor;
  }

  return 0;
}

int interaction::CUDA_Nonbonded_Set
::init(topology::Topology const & topo,
        configuration::Configuration const & conf,
        simulation::Simulation const & sim,
        std::ostream & os,
        bool quiet) {
  const int num_atoms = topo.num_atoms();

  m_storage.force.resize(num_atoms);
  m_storage_cuda.force.resize(num_atoms);
  m_longrange_storage.force.resize(num_atoms);
  m_longrange_storage_cuda.force.resize(num_atoms);
  m_storage.energies.
          resize(unsigned(conf.current().energies.bond_energy.size()),
          unsigned(conf.current().energies.kinetic_energy.size()));
  m_longrange_storage.energies.
          resize(unsigned(conf.current().energies.bond_energy.size()),
          unsigned(conf.current().energies.kinetic_energy.size()));

  m_longrange_storage_cuda.energies.
          resize(unsigned(conf.current().energies.bond_energy.size()),
          unsigned(conf.current().energies.kinetic_energy.size()));

  m_storage_cuda.energies.
          resize(unsigned(conf.current().energies.bond_energy.size()),
          unsigned(conf.current().energies.kinetic_energy.size()));

  m_storage.electric_field.resize(num_atoms);
  m_longrange_storage.electric_field.resize(num_atoms);

  DEBUG(10, "pairlist size: " << num_atoms);
  pairlist().resize(num_atoms);

  // check if we can guess the number of pairs
  const double vol = math::volume(conf.current().box, conf.boundary_type);
  if (vol) {
    const double c3 = sim.param().pairlist.cutoff_short *
            sim.param().pairlist.cutoff_short *
            sim.param().pairlist.cutoff_short;

    const unsigned int pairs =
            int(1.3 * num_atoms / vol * 4.0 / 3.0 * math::Pi * c3);

    if (!quiet)
      os << "\n\testimated pairlist size (per atom) : "
            << pairs << "\n";

    pairlist().reserve(pairs);
  }
  std::cout << "CUDA_Nonbonded_set:init" << std::endl;
  return 0;
}
