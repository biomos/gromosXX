/**
 * @file nonbonded_interaction.cc
 * template methods of Nonbonded_Interaction.
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

#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Nonbonded_Interaction::Nonbonded_Interaction(Pairlist_Algorithm *pa)
  : Interaction("NonBonded"),
    m_pairlist_algorithm(pa),
    m_parameter(),
    m_longrange_timing(0.0)
{
}

/**
 * Destructor.
 * @TODO change destruction of nonbonded set to be standard - conform!
 */
interaction::Nonbonded_Interaction::~Nonbonded_Interaction()
{
  DEBUG(4, "Nonbonded_Interaction::destructor");
  delete m_pairlist_algorithm;
  for(std::vector<Nonbonded_Set *>::iterator it = m_nonbonded_set.begin(), to = m_nonbonded_set.end();
      it != to;
      ++it)
    delete *it;
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::Nonbonded_Interaction::calculate_interactions(topology::Topology & topo,
							       configuration::Configuration & conf,
							       simulation::Simulation & sim)
{
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  const double nonbonded_start = util::now();

  // shared memory do this only once
  m_pairlist_algorithm->prepare(topo, conf, sim);

  std::vector<Nonbonded_Set *>::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();

#ifdef OMP
  int tid;
#pragma omp parallel private(tid)
  {
    tid = omp_get_thread_num();
    // calculate the corresponding interactions
    assert(m_nonbonded_set.size() > tid);
    m_nonbonded_set[tid]->calculate_interactions(topo, conf, sim,
						 tid, m_omp_num_threads);
  }
  
#else
  
  // have to do all from here (probably it's only one, could unite this,
  // but then maybe it's clearer like it is...
  
  for( ; it != to; ++it){
    (*it)->calculate_interactions(topo, conf, sim);
  }

#endif

  // add the forces, energies, virial...
  it = m_nonbonded_set.begin();

  const double ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & e = conf.current().energies;
  
  for( ; it != to; ++it){
    conf.current().force += (*it)->shortrange_storage().force;

    for(unsigned int i = 0; i < ljs; ++i){
      for(unsigned int j = 0; j < ljs; ++j){
      
	e.lj_energy[i][j] += 
	  (*it)->shortrange_storage().energies.lj_energy[i][j];
	e.crf_energy[i][j] += 
	  (*it)->shortrange_storage().energies.crf_energy[i][j];
      }
    }
    
    if (sim.param().pcouple.virial){
      DEBUG(7, "\tadd long range virial");

      for(unsigned int i=0; i<3; ++i){
	for(unsigned int j=0; j<3; ++j){

	  DEBUG(8, "set virial = " << (*it)->shortrange_storage().virial_tensor(i,j)
		<< "\tvirial = " << conf.current().virial_tensor(i,j));
	  
	  conf.current().virial_tensor(i,j) +=
	    (*it)->shortrange_storage().virial_tensor(i,j);
	}
      }
      
    }
  }
  
  if (sim.param().perturbation.perturbation){

    it = m_nonbonded_set.begin();
      
    for( ; it != to; ++it){
	
      configuration::Energy & pe = conf.current().perturbed_energy_derivatives;
	
      for(unsigned int i = 0; i < ljs; ++i){
	for(unsigned int j = 0; j < ljs; ++j){
      
	  pe.lj_energy[i][j] += 
	    (*it)->shortrange_storage().perturbed_energy_derivatives.lj_energy[i][j];
	  pe.crf_energy[i][j] += 
	    (*it)->shortrange_storage().perturbed_energy_derivatives.crf_energy[i][j];
	}
      }
    } // sets
  } // perturbation

  m_timing += util::now() - nonbonded_start;

  return 0;
  
}

/**
 * calculate the hessian for a given atom.
 */
int interaction::Nonbonded_Interaction::calculate_hessian(topology::Topology & topo,
							  configuration::Configuration & conf,
							  simulation::Simulation & sim,
							  unsigned int atom_i, unsigned int atom_j,
							  math::Matrix & hessian)
{
  std::vector<Nonbonded_Set *>::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();

  hessian = 0.0;
  math::Matrix h;

  for( ; it != to; ++it){
    (*it)->calculate_hessian(topo, conf, sim, atom_i, atom_j, h);

    for(unsigned int d1=0; d1 < 3; ++d1)
      for(unsigned int d2=0; d2 < 3; ++d2)
	hessian(d1,d2) += h(d1,d2);
  }
  return 0;
}

/**
 * initialize the arrays
 */
int interaction::Nonbonded_Interaction::init(topology::Topology const & topo,
					     configuration::Configuration const & conf,
					     simulation::Simulation const & sim,
					     bool quiet)
{

#ifdef OMP
  int tid;
#pragma omp parallel private(tid)
    {
      tid = omp_get_thread_num();
      if (tid == 0){
	m_omp_num_threads = omp_get_num_threads();
      }
    }
#else
    m_omp_num_threads = 1;
#endif

  DEBUG(15, "nonbonded_interaction::initialize");

  if (sim.param().perturbation.perturbation){
    m_nonbonded_set.
      resize(m_omp_num_threads, new Perturbed_Nonbonded_Set(*m_pairlist_algorithm, m_parameter));
  }
  else{
    m_nonbonded_set.
      resize(m_omp_num_threads, new Nonbonded_Set(*m_pairlist_algorithm, m_parameter));
  }
  
  std::vector<Nonbonded_Set *>::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();
  
  for( ; it != to; ++it){
    (*it)->init(topo, conf, sim, quiet);
  }

  return 0;
}


//***************************************************************************
// helper functions 
//***************************************************************************

void interaction::Nonbonded_Interaction::print_timing(std::ostream & os)
{
  os << "        "
     << std::setw(36) << std::left << name
     << std::setw(20) << m_timing << "\n"
     << "            "
     << std::setw(32) << std::left << "shortrange"
     << std::setw(20) 
     << m_timing - m_pairlist_algorithm->timing()
     << "\n"
     << "            "
     << std::setw(32) << std::left << "longrange"
    // << std::setw(20) << m_longrange_timing << "\n"
     << std::setw(20) << "not measured: too expensive" << "\n"
    
     << "            "
     << std::setw(32) << std::left << "pairlist"
     << std::setw(20) 
     << m_pairlist_algorithm->timing() - m_longrange_timing<< "\n"
     << "\n";
}

