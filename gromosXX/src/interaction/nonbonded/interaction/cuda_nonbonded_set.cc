/**
 * @file cuda_nonbonded_set.cc
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction_types.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../../../interaction/nonbonded/interaction/cuda_nonbonded_set.h"


#include "../../../interaction/nonbonded/interaction/latticesum.h"

#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"

#include "../../../math/periodicity.h"
#include "../../../math/volume.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_set.h"

#include "../../../util/debug.h"

#include "../../../configuration/energy.h"

#include "storage.h"
#include "../../../util/cycle_thread.h"

#ifdef HAVE_LIBCUDART
#include <cudaKernel.h>
#endif


#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE cuda


/**
 * Constructor.
 */
interaction::CUDA_Nonbonded_Set
::CUDA_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param,
		int rank, int num_threads,
        unsigned int gpu_id)
  : Nonbonded_Set(pairlist_alg, param, rank, num_threads),
    util::CycleThread(){

  // these values should be fine for a cutoff of 0.8 / 1.4. They have
  // to be increased for larger cutoffs. 
  estNeigh_long = 500;
  estNeigh_short = 400;

  // Which gpu for this set
  mygpu_id = gpu_id;
  
  // For m_param
  m_parameter = & param;

  DEBUG (8, "CUDA_Nonbonded_Set::constructor")
  //
}

/**
 * Destructor.
 */
interaction::CUDA_Nonbonded_Set::~CUDA_Nonbonded_Set(){
  DEBUG (8, "CUDA_Nonbonded_Set::destructor")
  terminate_cycle();
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::CUDA_Nonbonded_Set
::calculate_interactions(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(8, "CUDA_Nonbonded_Set::calculate_interactions");
  DEBUG(15, "Start Cycle for GPU: " << mygpu_id);
  if (mygpu_id == 0)
    m_pairlist_alg.timer().start("cuda set");
  do_cycle();
  if (mygpu_id == 0)
    m_pairlist_alg.timer().stop("cuda set");
  if (error) return 1;
  DEBUG(15, "Finished Cycle for GPU: " << mygpu_id);

  return 0;
}

/**
 * Initialize the CUDA nonbonded set
 */
int interaction::CUDA_Nonbonded_Set
::init(topology::Topology const & topo,
        configuration::Configuration const & conf,
        simulation::Simulation const & sim,
        std::ostream & os,
        bool quiet) {

  // Set the pointers and variables
  mytopo = (topology::Topology *) &topo;
  myconf = (configuration::Configuration *)&conf;
  mysim = (simulation::Simulation *) & sim;
  mystream = &os;
  amIquiet = quiet;

  const int num_atoms = topo.num_atoms();

  m_storage.force.resize(num_atoms);
  m_longrange_storage.force.resize(num_atoms);
  m_storage.energies.
          resize(unsigned(conf.current().energies.bond_energy.size()),
          unsigned(conf.current().energies.kinetic_energy.size()));
  m_longrange_storage.energies.
          resize(unsigned(conf.current().energies.bond_energy.size()),
          unsigned(conf.current().energies.kinetic_energy.size()));

  start();
  pthread_barrier_wait(&barrier_init);
  return 0;
}

/**
 * calculate constants needed for the main calculations
 */
void interaction::CUDA_Nonbonded_Set::init_run() {

#ifdef HAVE_LIBCUDART
  /////////////////////////////////
  // Calculate constants for the solvent
  //
  DEBUG(8, "CUDA_Nonbonded_Set::run")

          unsigned int nAtomsPerSolvMol = mytopo->num_solvent_atoms() / mytopo->num_solvent_molecules(0);
  cudakernel::lj_crf_parameter * pLj_crf;
  pLj_crf = (cudakernel::lj_crf_parameter *)malloc(nAtomsPerSolvMol * nAtomsPerSolvMol * sizeof (cudakernel::lj_crf_parameter));

  unsigned int solvIdx = mytopo->num_solute_atoms();
  for (unsigned int i = 0; i < nAtomsPerSolvMol; i++) {
    for (unsigned int j = 0; j < nAtomsPerSolvMol; j++) {
      const lj_parameter_struct & lj = m_parameter->lj_parameter(mytopo->iac(solvIdx + i), mytopo->iac(solvIdx + j));
      pLj_crf[i * nAtomsPerSolvMol + j].c12 = lj.c12;
      pLj_crf[i * nAtomsPerSolvMol + j].c6 = lj.c6;
      pLj_crf[i * nAtomsPerSolvMol + j].q = math::four_pi_eps_i * mytopo->charge(solvIdx + i) * mytopo->charge(solvIdx + j);

    }
  }

  double m_cut3i, m_crf, m_crf_cut, m_crf_cut3i, m_crf_2cut3i;
  double rf_cutoff = mysim->param().nonbonded.rf_cutoff;
  double epsilon = mysim->param().nonbonded.epsilon;
  double rf_epsilon = mysim->param().nonbonded.rf_epsilon;
  double rf_kappa = mysim->param().nonbonded.rf_kappa;

  m_cut3i = 1.0 / (rf_cutoff * rf_cutoff * rf_cutoff);
  m_crf = 2 * (epsilon - rf_epsilon)*(1.0 + rf_kappa * rf_cutoff) - rf_epsilon *
          (rf_kappa * rf_cutoff * rf_kappa * rf_cutoff);
  m_crf /= (epsilon + 2 * rf_epsilon)*(1.0 + rf_kappa * rf_cutoff) +
          rf_epsilon * (rf_kappa * rf_cutoff * rf_kappa * rf_cutoff);
  m_crf_cut3i = m_crf*m_cut3i;

  m_crf_2cut3i = m_crf_cut3i / 2.0;
  m_crf_cut = (1 - m_crf / 2.0) / rf_cutoff;

  // check pressure coupling
  // anisotropic is ok now, because we can now handle different box edges
  // (only for a rectangular box; below we exit if the box is not rectangular)
  // so, only complain here if the pressure coupling is fully-anisotropic
  if (mysim->param().pcouple.scale != math::pcouple_off &&
          mysim->param().pcouple.scale == math::pcouple_full_anisotropic) {
//          mysim->param().pcouple.scale != math::pcouple_isotropic) {
//    io::messages.add("CUDA solvent doesn't support anisotropic pressure scaling.",
    io::messages.add("CUDA solvent doesn't support fully-anisotropic pressure scaling.",
            "CUDA_Nonbonded", io::message::error);
    return;
  }

  if (myconf->boundary_type !=  math::rectangular) { 
      DEBUG(9, "BOX IS NOT RECTANGULAR!!!")
      io::messages.add("Box is not rectangular!",
            "CUDA_Nonbonded", io::message::error);
    return;
  }

  //
  // calculation end here
  /////////////////////////////////////

  DEBUG(9, "CUDA_Nonbonded_Set::run : cudaInit")
  DEBUG(11, "mygpu_id: " << mygpu_id << " of " << mysim->param().innerloop.number_gpus << ". device number: " << mysim->param().innerloop.gpu_device_number.at(mygpu_id))
  gpu_stat = cudakernel::cudaInit
          (
          //mygpu_id, /*sim.param().innerloop.cuda_device,*/
          mysim->param().innerloop.gpu_device_number.at(mygpu_id),
          mytopo->num_solvent_atoms(),
          mysim->param().pairlist.cutoff_short,
          mysim->param().pairlist.cutoff_long,
          myconf->current().box(0)(0),
          myconf->current().box(1)(1),
          myconf->current().box(2)(2),
          mytopo->num_solvent_atoms() / mytopo->num_solvent_molecules(0),
          //cuda_nbs->estNeigh_short,
          //cuda_nbs->estNeigh_long,
          this->estNeigh_short,
          this->estNeigh_long,
          m_crf_2cut3i,
          m_crf_cut,
          m_crf_cut3i,
          pLj_crf,
          mysim->param().innerloop.number_gpus,
          mygpu_id,
          &error
          );
  if (error) {
    std::ostringstream msg;
    msg << "Cannot initialize nonbonded interaction on GPU " << mysim->param().innerloop.gpu_device_number.at(mygpu_id);
    io::messages.add(msg.str(), "CUDA_Nonbonded_Set", io::message::error);
    return;
  }
  free(pLj_crf);
  return;
#endif
}

// ----------------------------------------------------------------------------
// CYCLE

/**
 * what the thread does: calculate the solvent forces, energies and virals
 */
void interaction::CUDA_Nonbonded_Set::cycle() {
  DEBUG(6, "CUDA_Nonbonded_Set::cycle start calculations. GPU: " << mygpu_id);
  // this whole function is remove is the lib is not present.

#ifdef HAVE_LIBCUDART

  m_storage.zero();
  const bool pairlist_update = !(mysim->steps() % mysim->param().pairlist.skip_step);
  if (mygpu_id == 0)  
    m_pairlist_alg.timer().start("GPU data copy");
  
  error += cudakernel::cudaCopyPositions(&myconf->current().pos(mytopo->num_solute_atoms())(0), gpu_stat);
  DEBUG(15, "myconf->current().pos(mytopo->num_solute_atoms())(0) = " << myconf->current().pos(mytopo->num_solute_atoms())(0));
  

  // copy the box if pressure is coupled
  if (mysim->param().pcouple.scale != math::pcouple_off) {
    error += cudakernel::cudaCopyBox(gpu_stat, myconf->current().box(0)(0), myconf->current().box(1)(1), myconf->current().box(2)(2));
  }
  if (mygpu_id == 0)  
    m_pairlist_alg.timer().stop("GPU data copy");
  

  if (pairlist_update) {
    DEBUG(6, "\tdoing longrange...");

    m_longrange_storage.zero();
    if (mygpu_id == 0)
      m_pairlist_alg.timer().start("pairlist cuda");
    cudakernel::cudaCalcPairlist(gpu_stat);
    if (mygpu_id == 0)
      m_pairlist_alg.timer().stop("pairlist cuda");
  }

  if (pairlist_update) {
    const int egroup = mytopo->atom_energy_group(mytopo->num_solute_atoms());
    DEBUG(15, "egroup = " << egroup);
    DEBUG(8, "doing longrange calculation");

    // in case of LS only LJ are calculated

    if (mygpu_id == 0)
      m_pairlist_alg.timer().start("longrange-cuda");
    double * For = &m_longrange_storage.force(mytopo->num_solute_atoms())(0);
    DEBUG(15, "mytopo->num_solute_atoms() = " << mytopo->num_solute_atoms());
    double * Vir = &m_longrange_storage.virial_tensor(0, 0);
    double * e_lj = &m_longrange_storage.energies.lj_energy[egroup][egroup];
    double * e_crf = &m_longrange_storage.energies.crf_energy[egroup][egroup];
    DEBUG(15, "For = " << *For);
    DEBUG(15, "Vir = " << *Vir);
    DEBUG(15, "e_lj = " << *e_lj);
    DEBUG(15, "e_crf = " << *e_crf);
    error += cudakernel::cudaCalcForces(For, Vir, e_lj, e_crf, true, gpu_stat);
    if (mygpu_id == 0)
      m_pairlist_alg.timer().stop("longrange-cuda");
  }
  // calculate forces / energies
  DEBUG(6, "\tshort range interactions");

  if (mygpu_id == 0)
    m_pairlist_alg.timer().start("shortrange-cuda");
  double * For = &m_storage.force(mytopo->num_solute_atoms())(0);
  double * Vir = &m_storage.virial_tensor(0, 0);
  const int egroup = mytopo->atom_energy_group(mytopo->num_solute_atoms());
  double * e_lj = &m_storage.energies.lj_energy[egroup][egroup];
  double * e_crf = &m_storage.energies.crf_energy[egroup][egroup];
  error += cudakernel::cudaCalcForces(For, Vir, e_lj, e_crf, false, gpu_stat);
  if (mygpu_id == 0)
    m_pairlist_alg.timer().stop("shortrange-cuda");
  if (m_rank == 0 && error) {
    io::messages.add("GPU: cannot calculate forces", io::message::critical);
    return;
  }

  // add long-range force
  DEBUG(6, "\t(set) add long range forces");

  m_storage.force += m_longrange_storage.force;

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

  // add longrange virial
  if (mysim->param().pcouple.virial) {
    DEBUG(6, "\t(set) add long range virial");

    m_storage.virial_tensor += m_longrange_storage.virial_tensor;
  }
#endif
}

// ----------------------------------------------------------------------------
// CLEANUP

/**
 * Clean up after cycle
 */
void interaction::CUDA_Nonbonded_Set::end_run() {
  DEBUG(15, "CUDA_Nonbonded_Set: Cleaning up for GPU: " << mygpu_id);
#ifdef HAVE_LIBCUDART
  if (cudakernel::CleanUp(gpu_stat))
    io::messages.add("GPU cleanup failed", io::message::critical);
#endif
}
