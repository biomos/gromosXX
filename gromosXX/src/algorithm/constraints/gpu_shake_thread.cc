/**
 * @file gpu_shake_thread.cc
 * Source file, which implements pthread for M_SHAKE on the GPU
 * Similar to cudaNonbondedSet
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../math/periodicity.h"

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

#include "../../math/gmath.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * Constructor
 */
algorithm::GPU_Shake_Thread::GPU_Shake_Thread(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        math::GenericMatrix<double> factor,
        math::Vec mass,
        math::Vec constr_length2,
        double tolerance,
        unsigned int num_gpus,
        unsigned int gpu_id,
        std::ostream & os,
        bool quiet) : util::CycleThread(),
                               mytopo(&topo), myconf(&conf), mysim(&sim),
                               factor(factor), mass(mass), constr_length2(constr_length2),
                               m_tolerance(tolerance),
                               num_gpus(num_gpus), gpu_id(gpu_id),
                               mystream(&os), amIquiet(quiet) {
  shake_fail_mol = -1;
  DEBUG(5, "Constructor of GPU_Shake_Thread " << gpu_id);
  start();
  pthread_barrier_wait(&barrier_init);
}

/**
 * Destructor
 */
algorithm::GPU_Shake_Thread::~GPU_Shake_Thread() {
  DEBUG(5, "Destructor of GPU_Shake_Thread " << gpu_id)
  terminate_cycle();
}


/**
 * Apply the algorithm
 */
int algorithm::GPU_Shake_Thread::apply(){
  DEBUG(10, "GPU_Shake_Thread : Do Cycle (" << gpu_id << ")")
  do_cycle();

  return shake_fail_mol;
}


/**
 * Initialize the GPU
 */
void algorithm::GPU_Shake_Thread::init_run() {
  DEBUG(10, "M_SHAKE : Init : Number of GPUs: " << mysim->param().constraint.solvent.number_gpus << " ID: " << mysim->param().constraint.solvent.gpu_device_number.at(gpu_id));
#ifdef HAVE_LIBCUDART
  int dev = mysim->param().constraint.solvent.gpu_device_number[gpu_id];
  int error = 0;

  if (mysim->param().innerloop.method == simulation::sla_cuda && dev == -1) {
    if (gpu_id < mysim->param().innerloop.number_gpus) {
      dev = mysim->param().innerloop.gpu_device_number[gpu_id];
    }
  }
  gpu_stat = cudakernel::cudaInitGPU_Shake(
          dev,
          &constr_length2(0),
          &factor(0, 0), &mass(0),
          m_tolerance,
          mysim->param().constraint.solvent.number_gpus,
          gpu_id,
          mytopo->num_solvent_atoms(0),
          mytopo->num_solvent_molecules(0),
          &error);
  mysim->param().constraint.solvent.gpu_device_number[gpu_id] = dev;
  if (error) {
    std::ostringstream msg;
    msg << "Cannot initialize SHAKE on GPU " << dev;
    io::messages.add(msg.str(), "GPU_Shake_Thread", io::message::error);
    return;
  }
#endif
  return;
}

/**
 * What the cycle does: calculating the constraints
 */
void algorithm::GPU_Shake_Thread::cycle(){

    shake_fail_mol = -1;

    DEBUG(10, "GPU_Shake_Thread : Cycle : Calculate Constraints")
#ifdef HAVE_LIBCUDART
    cudakernel::cudaGPU_Shake(&(myconf->current().pos(mytopo->num_solute_atoms())(0)),
          &(myconf->old().pos(mytopo->num_solute_atoms())(0)),
          shake_fail_mol, gpu_stat);
#endif
}


/**
 * Clean up
 */
void algorithm::GPU_Shake_Thread::end_run() {
  DEBUG(10, "GPU_Shake_Thread : Clean up");
#ifdef HAVE_LIBCUDART
  cudakernel::CleanUp(gpu_stat);
#endif
}


