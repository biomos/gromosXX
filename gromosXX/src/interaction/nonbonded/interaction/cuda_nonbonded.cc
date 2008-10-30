/**
 * @file cuda_nonbonded.cc
 * methods of CUDA_Nonbonded
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <simulation/parameter.h>
#include <configuration/configuration.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>

#include <interaction/nonbonded/interaction/storage.h>

#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>

#include <interaction/nonbonded/interaction/nonbonded_set.h>
#include <interaction/nonbonded/interaction/cuda_nonbonded_set.h>
#include <interaction/nonbonded/interaction/cuda_nonbonded.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/nonbonded_interaction.h>


#include <util/debug.h>
#include <util/error.h>

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
interaction::CUDA_Nonbonded::CUDA_Nonbonded(Pairlist_Algorithm *pa)
: Nonbonded_Interaction(pa)
{}

/**
 * Destructor.
 */
interaction::CUDA_Nonbonded::~CUDA_Nonbonded() 
{
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::CUDA_Nonbonded::
calculate_interactions(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) 
{
  DEBUG(4, "CUDA_NonBonded::calculate_interactions");
  
  m_timer.start();
  m_pairlist_algorithm->prepare(topo, conf, sim);
 
  if(m_nonbonded_set[0]->calculate_interactions(topo, conf, sim) > 0)
	return 1;

  store_set_data(topo,conf,sim);
  
  DEBUG(6, "CUDA_NonBonded_Interaction::calculate_interactions done");
  m_timer.stop();
  
  return 0;
  
}

int interaction::CUDA_Nonbonded::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) 
{
#ifdef HAVE_LIBCUKERNEL
  m_pairlist_algorithm->init(topo, conf, sim, os, quiet);
  CUDA_Nonbonded_Set * cuda_nbs = new CUDA_Nonbonded_Set(*m_pairlist_algorithm,m_parameter,0,1);

  unsigned int nAtomsPerSolvMol=topo.num_solvent_atoms()/topo.num_solvent_molecules(0);
  cudakernel::lj_crf_parameter * pLj_crf;
  pLj_crf=(cudakernel::lj_crf_parameter *)malloc(nAtomsPerSolvMol*nAtomsPerSolvMol*sizeof(cudakernel::lj_crf_parameter));
  
  unsigned int solvIdx=topo.num_solute_atoms();
  for(unsigned int i=0;i<nAtomsPerSolvMol;i++){
            for(unsigned int j=0;j<nAtomsPerSolvMol;j++){
                    const lj_parameter_struct & lj=m_parameter.lj_parameter(topo.iac(solvIdx+i),topo.iac(solvIdx+j));
                    pLj_crf[i*nAtomsPerSolvMol+j].c12 = lj.c12;
                    pLj_crf[i*nAtomsPerSolvMol+j].c6 = lj.c6;
                    pLj_crf[i*nAtomsPerSolvMol+j].q = topo.charge(solvIdx+i)*topo.charge(solvIdx+j);
                    
            }
  }
  
  double m_cut3i, m_crf, m_crf_cut,m_crf_cut3i,m_crf_2cut3i;
  double rf_cutoff=sim.param().nonbonded.rf_cutoff;
  double epsilon=sim.param().nonbonded.epsilon;
  double rf_epsilon=sim.param().nonbonded.rf_epsilon;
  double rf_kappa=sim.param().nonbonded.rf_kappa;
  
  m_cut3i=1.0/(rf_cutoff*rf_cutoff*rf_cutoff);
  m_crf=2*(epsilon-rf_epsilon)*(1.0+rf_kappa*rf_cutoff)-rf_epsilon*
		(rf_kappa*rf_cutoff*rf_kappa*rf_cutoff);
  m_crf /= (epsilon+2*rf_epsilon)*(1.0+rf_kappa*rf_cutoff)+
		rf_epsilon*(rf_kappa*rf_cutoff*rf_kappa*rf_cutoff);
  m_crf_cut3i=m_crf*m_cut3i;

  m_crf_2cut3i=m_crf_cut3i / 2.0;
  m_crf_cut=(1-m_crf/2.0)/rf_cutoff;

  cudakernel::cudaInit
	(
		topo.num_solvent_atoms(),
		sim.param().pairlist.cutoff_short,
		sim.param().pairlist.cutoff_long,
		conf.current().box(0)(0),
		topo.num_solvent_atoms()/topo.num_solvent_molecules(0),
		topo.num_solute_atoms(),
                cuda_nbs->estNeigh_long,
                cuda_nbs->estNeigh_short,
		m_crf_2cut3i,
		m_crf_cut,
		m_crf_cut3i,
                pLj_crf
	);

  //} initialize CUDA
  free(pLj_crf);
  DEBUG(15, "CUDA_NonBonded::initialize");
  m_nonbonded_set.clear();
  
  m_nonbonded_set.push_back(cuda_nbs);
  
  m_nonbonded_set[0]->init(topo, conf, sim, os, quiet);
#else
  io::messages.add("CUDA kernel initialized but no CUDA library.",
          "CUDA_Nonbonded", io::message::critical);
#endif

  return 0;

  
}
