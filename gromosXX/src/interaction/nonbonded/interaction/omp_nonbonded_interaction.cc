/**
 * @file omp_nonbonded_interaction.cc
 * methods of OMP_Nonbonded_Interaction.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>

#include <interaction/nonbonded/interaction/storage.h>

#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>
#include <interaction/nonbonded/interaction/nonbonded_set.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_term.h>

#include <interaction/nonbonded/interaction/perturbed_nonbonded_pair.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_set.h>

#include <interaction/nonbonded/interaction/nonbonded_interaction.h>
#include <interaction/nonbonded/interaction/omp_nonbonded_interaction.h>

#include <util/debug.h>
#include <util/error.h>

#ifdef OMP
#include <omp.h>
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

  assert((sim.param().force.spc_loop <= 0) || 
	 (!sim.param().pairlist.grid && !sim.param().pairlist.atomic_cutoff));

  const double nonbonded_start = util::now();

  // check if we want to calculate nonbonded
  // might not be necessary if multiple time-stepping is enabled

  int steps = sim.param().multistep.steps;
  if (steps == 0) steps = 1;
  
  if ((sim.steps() % steps) == 0){
    // std::cout << "MULTISTEP: full calculation\n";

    ////////////////////////////////////////////////////
    // multiple unit cell
    ////////////////////////////////////////////////////
    if (sim.param().multicell.multicell){
      
      DEBUG(6, "nonbonded: MULTICELL");
      configuration::Configuration exp_conf;
      expand_configuration(topo, conf, sim, exp_conf);
      DEBUG(7, "\tmulticell conf: pos.size()=" << exp_conf.current().pos.size());
      
      // shared memory do this only once
      m_pairlist_algorithm->prepare(topo.multicell_topo(), exp_conf, sim);
      
#ifdef OMP
      int tid;
#pragma omp parallel private(tid)
      {
	tid = omp_get_thread_num();
	// calculate the corresponding interactions
	assert(m_nonbonded_set.size() > tid);
	DEBUG(8, "calculating nonbonded interactions (thread " 
	      << tid << " of " << m_set_size << ")");
	
	m_nonbonded_set[tid]->calculate_interactions(topo.multicell_topo(), 
						     exp_conf, sim);
      }
      
#else
    
      std::cerr << "using OMP code without OMP defined..." << std::endl;
      return E_ILLEGAL;
      
#endif
    }
    else{ // no MULTICELL
    
      // shared memory do this only once
      m_pairlist_algorithm->prepare(topo, conf, sim);
      
#ifdef OMP
      int tid;
#pragma omp parallel private(tid)
      {
	tid = omp_get_thread_num();
	// calculate the corresponding interactions
	assert(m_nonbonded_set.size() > tid);
	DEBUG(8, "calculating nonbonded interactions (thread " 
	      << tid << " of " << m_set_size << ")");
	
	m_nonbonded_set[tid]->calculate_interactions(topo, conf, sim);
      }
      
#else
      std::cerr << "using OMP code without OMP defined..." << std::endl;
      return E_ILLEGAL;
#endif
    }
    
    ////////////////////////////////////////////////////
    // end of multiple time stepping: calculate
    ////////////////////////////////////////////////////
  }
  else{
    // std::cout << "MULTISTEP: no recalculation...\n";
  }
  
  
  DEBUG(6, "sets are done, adding things up...");
  store_set_data(topo, conf, sim);

  DEBUG(6, "Nonbonded_Interaction::calculate_interactions done");

  m_timing += util::now() - nonbonded_start;

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
  int tid;
#pragma omp parallel private(tid)
    {
      tid = omp_get_thread_num();
      if (tid == 0){
	m_set_size = omp_get_num_threads();
      }
    }
#else
    std::cerr << "OMP not defined, why are we here???" << std::endl;
    return E_ILLEGAL;
#endif
    
    return Nonbonded_Interaction::init(topo, conf, sim, os, quiet);
}


//***************************************************************************
// helper functions 
//***************************************************************************


