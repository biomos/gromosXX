/**
 * @file mpi_nonbonded_master.cc
 * methods of MPI_Nonbonded_Master
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
#include <interaction/nonbonded/interaction/mpi_nonbonded_master.h>

#include <util/debug.h>
#include <util/error.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::MPI_Nonbonded_Master::MPI_Nonbonded_Master(Pairlist_Algorithm *pa)
  : Nonbonded_Interaction(pa)
{
}

/**
 * Destructor.
 * @bug change destruction of nonbonded set to be standard - conform!
 */
interaction::MPI_Nonbonded_Master::~MPI_Nonbonded_Master()
{
  DEBUG(7, "MPI_Nonbonded_Master::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::MPI_Nonbonded_Master::
calculate_interactions(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim)
{
  DEBUG(4, "MPI_Nonbonded_Master::calculate_interactions");

  assert((!sim.param().force.spc_loop) || 
	 (!sim.param().pairlist.grid && !sim.param().pairlist.atomic_cutoff));

  const double nonbonded_start = util::now();

  // do this on the master and on the slaves...
  m_pairlist_algorithm->prepare(topo, conf, sim);
  
#ifdef XXMPI
  int rank = MPI::COMM_WORLD.Get_rank();
  int num_threads = MPI::COMM_WORLD.Get_size();

  // distribute the positions

  // if this is true, the positions (and forces) can be copied
  // directly from theVArray
  assert(((double *) &conf.current().pos(conf.current().pos.size()-1)(0)) -
	 ((double *) &conf.current().pos(0)(0)) == (conf.current().pos.size() - 1)*3);

  MPI::COMM_WORLD.Bcast(&conf.current().pos(0)(0),
			conf.current().pos.size() * 3, 
			MPI::DOUBLE,
			0);

  // calculate interactions for our rank
  DEBUG(8, "calculating nonbonded interactions (thread " 
	<< rank << " of " << num_threads << ")");
    
  m_nonbonded_set[0]->calculate_interactions(topo, conf, sim);

  // collect the forces, energies, energy-derivatives, virial
  // MPI::IN_PLACE ???
  m_storage.force = m_nonbonded_set[0]->shortrange_storage().force;
  MPI::COMM_WORLD.Reduce(&m_storage.force(0)(0),
			 &m_nonbonded_set[0]->shortrange_storage().force(0),
			 m_nonbonded_set[0]->shortrange_storage().force.size() * 3,
			 MPI::DOUBLE,
			 MPI::SUM,
			 0);

  const unsigned int ljs = conf.current().energies.lj_energy.size();
  std::vector<double> lj_scratch(ljs*ljs), rlj_scratch(ljs*ljs);
  std::vector<double> crf_scratch(ljs*ljs), rcrf_scratch(ljs*ljs);
  for(unsigned int i = 0; i < ljs; ++i){
    for(unsigned int j = 0; j < ljs; ++j){
      lj_scratch[i*ljs + j] = 
	m_nonbonded_set[0]->shortrange_storage().energies.lj_energy[i][j];
      crf_scratch[i*ljs + j] = 
	m_nonbonded_set[0]->shortrange_storage().energies.crf_energy[i][j];
    }
  }
  MPI::COMM_WORLD.Reduce(&lj_scratch[0],
			 &rlj_scratch[0],
			 ljs * ljs,
			 MPI::DOUBLE,
			 MPI::SUM,
			 0);
  MPI::COMM_WORLD.Reduce(&crf_scratch[0],
			 &rcrf_scratch[0],
			 ljs * ljs,
			 MPI::DOUBLE,
			 MPI::SUM,
			 0);
  
  for(unsigned int i = 0; i < ljs; ++i){
    for(unsigned int j = 0; j < ljs; ++j){
      m_nonbonded_set[0]->shortrange_storage().energies.lj_energy[i][j] = rlj_scratch[i*ljs + j];
      m_nonbonded_set[0]->shortrange_storage().energies.crf_energy[i][j] = rcrf_scratch[i*ljs + j];
    }
  }

  if (sim.param().pcouple.virial){

    m_storage.virial_tensor = m_nonbonded_set[0]->shortrange_storage().virial_tensor;
    double * dvt = &m_storage.virial_tensor(0,0);
    double * dvt2 = &m_nonbonded_set[0]->shortrange_storage().virial_tensor(0,0);
    MPI::COMM_WORLD.Reduce(dvt,
			   dvt2,
			   9,
			   MPI::DOUBLE,
			   MPI::SUM,
			   0);
  }
  
  if (sim.param().perturbation.perturbation){

    for(unsigned int i = 0; i < ljs; ++i){
      for(unsigned int j = 0; j < ljs; ++j){
	lj_scratch[i*ljs + j] = 
	  m_nonbonded_set[0]->shortrange_storage().perturbed_energy_derivatives.lj_energy[i][j];
	crf_scratch[i*ljs + j] = 
	  m_nonbonded_set[0]->shortrange_storage().perturbed_energy_derivatives.crf_energy[i][j];
      }
    }
    MPI::COMM_WORLD.Reduce(&lj_scratch[0],
			   &rlj_scratch[0],
			   ljs * ljs,
			   MPI::DOUBLE,
			   MPI::SUM,
			   0);
    MPI::COMM_WORLD.Reduce(&crf_scratch[0],
			   &rcrf_scratch[0],
			   ljs * ljs,
			   MPI::DOUBLE,
			   MPI::SUM,
			   0);
  
    for(unsigned int i = 0; i < ljs; ++i){
      for(unsigned int j = 0; j < ljs; ++j){
	m_nonbonded_set[0]->shortrange_storage().perturbed_energy_derivatives.lj_energy[i][j]
	  = rlj_scratch[i*ljs + j];
	m_nonbonded_set[0]->shortrange_storage().perturbed_energy_derivatives.crf_energy[i][j]
	  = rcrf_scratch[i*ljs + j];
      }
    }
  }
  
#else
  std::cerr << "using MPI code without MPI defined..." << std::endl;
  return E_ILLEGAL;
#endif
  
  DEBUG(6, "sets are done, adding things up...");
  store_set_data(topo, conf, sim);
  
  DEBUG(6, "Nonbonded_Interaction::calculate_interactions done");
  
  m_timing += util::now() - nonbonded_start;
  
  return 0;
  
}

/**
 * initialize the arrays
 * need to override to pass MPI rank, size to nonbonded set
 */
int interaction::MPI_Nonbonded_Master::init(topology::Topology & topo,
					    configuration::Configuration & conf,
					    simulation::Simulation & sim,
					    std::ostream & os,
					    bool quiet)
{
#ifdef XXMPI
  int rank = MPI::COMM_WORLD.Get_rank();
  int num_threads = MPI::COMM_WORLD.Get_size();

  // initialise the pairlist...
  m_pairlist_algorithm->init(topo, conf, sim, os, quiet);

  DEBUG(15, "MPI_Nonbonded_Master::initialize");
  m_nonbonded_set.clear();

  if (sim.param().perturbation.perturbation){

    // only one set per MPI process
    m_nonbonded_set.push_back(new Perturbed_Nonbonded_Set(*m_pairlist_algorithm,
							  m_parameter, rank, num_threads));
  }
  else{
    // only one set per MPI process
    m_nonbonded_set.push_back(new Nonbonded_Set(*m_pairlist_algorithm, 
						m_parameter, rank, num_threads));
  }
  
  if (sim.param().multicell.multicell)
    m_nonbonded_set[0]->init(topo.multicell_topo(), conf, sim, os, quiet);
  else
    m_nonbonded_set[0]->init(topo, conf, sim, os, quiet);

  m_storage.force.resize(topo.num_atoms());
  m_storage.energies.
    resize(unsigned(conf.current().energies.bond_energy.size()),
	   unsigned(conf.current().energies.kinetic_energy.size()));

  check_spc_loop(topo, conf, sim, os, quiet);
  return 0;
  
#else
  std::cerr << "MPI: MPI_Nonbonded_Master::init but MPI not defined" << std::endl;
  return 1;
#endif

}


//***************************************************************************
// helper functions 
//***************************************************************************

