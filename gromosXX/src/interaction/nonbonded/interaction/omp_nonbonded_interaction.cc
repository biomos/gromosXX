/**
 * @file omp_nonbonded_interaction.cc
 * methods of OMP_Nonbonded_Interaction.
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction.h"
#include "../../../interaction/interaction_types.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"

#include "../../../interaction/nonbonded/interaction/storage.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_set.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_term.h"

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_pair.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_set.h"
#include "../../../interaction/nonbonded/interaction/cuda_nonbonded_set.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_interaction.h"
#include "../../../interaction/nonbonded/interaction/omp_nonbonded_interaction.h"

#include "../../../util/debug.h"
#include "../../../util/error.h"
#include "../../../util/thread.h"

#include "../../../math/boundary_checks.h"

#ifdef OMP
#include <omp.h>
#endif

#ifdef HAVE_LIBCUDART
#include <cudaKernel.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded


/**
 * Constructor.
 */
interaction::OMP_Nonbonded_Interaction::OMP_Nonbonded_Interaction(Pairlist_Algorithm *pa)
  : Nonbonded_Interaction(pa)
{
}

/**
 * Destructor.
 * @bug change destruction of nonbonded set to be standard - conform!
 */
interaction::OMP_Nonbonded_Interaction::~OMP_Nonbonded_Interaction()
{
  DEBUG(7, "OMP_Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::OMP_Nonbonded_Interaction::
calculate_interactions(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim)
{
  DEBUG(4, "OMP_Nonbonded_Interaction::calculate_interactions");

  m_timer.start(sim);

  // check if we want to calculate nonbonded
  // might not be necessary if multiple time-stepping is enabled

  int steps = sim.param().multistep.steps;
  if (steps == 0) steps = 1;

  configuration::Configuration *p_conf = &conf;
  topology::Topology *p_topo = &topo;

  if ((sim.steps() % steps) == 0) {
    // std::cout << "MULTISTEP: full calculation\n";

    ////////////////////////////////////////////////////
    // multiple unit cell
    ////////////////////////////////////////////////////
    if (sim.param().multicell.multicell) {
      DEBUG(6, "nonbonded: MULTICELL");
      p_conf = m_exp_conf;
      p_topo = &topo.multicell_topo();
      expand_configuration(topo, conf, sim, *p_conf);
      if (!math::boundary_check_cutoff(p_conf->current().box, p_conf->boundary_type,
              sim.param().pairlist.cutoff_long)) {
        io::messages.add("box is too small: not twice the cutoff!",
                "configuration", io::message::error);
        return 1;
      }
      DEBUG(6, "\tmulticell conf: pos.size()=" << p_conf->current().pos.size());
    }

    // shared memory do this only once
    m_pairlist_algorithm->prepare(*p_topo, *p_conf, sim);

    int error = 0;
    #ifdef OMP
      omp_set_num_threads(m_nonbonded_set.size());
      #pragma omp parallel reduction(+:error)
      {
        unsigned int tid = omp_get_thread_num();
        // calculate the corresponding interactions
        assert(m_nonbonded_set.size() > tid);
        DEBUG(4, "calculating nonbonded interactions (thread "
              << tid << " of " << m_set_size << ")");

        error += m_nonbonded_set[tid]->calculate_interactions(*p_topo, *p_conf, sim);

      }

    #else
      std::cerr << "using OMP code without OMP defined..." << std::endl;
      return E_ILLEGAL;
    #endif
    if (error) return 1;
    
    ////////////////////////////////////////////////////
    // end of multiple time stepping: calculate
    ////////////////////////////////////////////////////
  }
  else{
    // std::cout << "MULTISTEP: no recalculation...\n";
  }
  
  
  DEBUG(6, "sets are done, adding things up...");
  store_set_data(*p_topo, *p_conf, sim);

  if (sim.param().multicell.multicell) {
    reduce_configuration(topo, conf, sim, *p_conf);
  }

  ////////////////////////////////////////////////////
  // printing pairlist
  ////////////////////////////////////////////////////
  if (sim.param().pairlist.print &&
      (!(sim.steps() % sim.param().pairlist.skip_step))) {
    DEBUG(7, "print pairlist...");
    std::cerr << "printing pairlist!" << std::endl;
    print_pairlist(*p_topo, *p_conf, sim);
  }
  
  DEBUG(6, "Nonbonded_Interaction::calculate_interactions done");

  m_timer.stop();

  return 0;
  
}

/**
 * initialize the arrays
 */
int interaction::OMP_Nonbonded_Interaction::init(topology::Topology & topo,
						 configuration::Configuration & conf,
						 simulation::Simulation & sim,
						 std::ostream & os,
						 bool quiet)
{

 

  // OpenMP parallelization
#ifdef OMP
  unsigned int number_of_cpus = 0;
  int result = 0;
  
  #pragma omp parallel
  {
    number_of_cpus = omp_get_num_threads();
    unsigned int tid = omp_get_thread_num();
    if (tid == 0)
      m_set_size = number_of_cpus;
  }
  result += Nonbonded_Interaction::init(topo, conf, sim, os, quiet);

  // Increase the number of threads to include the GPUs if CUDA enabled
  if (sim.param().innerloop.method == simulation::sla_cuda && result == 0) {
    omp_set_num_threads(number_of_cpus + sim.param().innerloop.number_gpus);

    #pragma omp parallel
    {
      unsigned int tid = omp_get_thread_num();
      DEBUG(10, "tid: " << tid);
      // For the GPUs
      if (tid >= number_of_cpus) {
        unsigned int gpu_tid = tid - number_of_cpus;
        DEBUG(9, "OMP: CUDA_set, gpu_tid: " << gpu_tid);

        CUDA_Nonbonded_Set * cuda_nbs = new CUDA_Nonbonded_Set(*m_pairlist_algorithm, m_parameter, 0, 1, gpu_tid);

        cuda_nbs->init(topo, conf, sim, os, quiet);

        DEBUG(10, "CUDA_NonBonded::initialize");

        #pragma omp critical
        {
          m_nonbonded_set.push_back(cuda_nbs);
        }
      }
    } // CUDA enabled
  }
  DEBUG(10, "OMP: result: " << result);

  return result;
#else
  std::cerr << "OMP not defined, why are we here???" << std::endl;
  return E_ILLEGAL;
#endif
}


