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
  DEBUG(7, "Nonbonded_Interaction::destructor");
  delete m_pairlist_algorithm;
  DEBUG(12, "pairlist algorithm destroyed");

  for(unsigned int i=0; i < m_nonbonded_set.size(); ++i){
    DEBUG(12, "deleting set " << i);
    delete m_nonbonded_set[i];
    m_nonbonded_set[i] = NULL;
  }
  
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::Nonbonded_Interaction::
calculate_interactions(topology::Topology & topo,
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
    DEBUG(8, "calculating nonbonded interactions (thread " 
	  << tid << " of " << m_omp_num_threads << ")");
    
    m_nonbonded_set[tid]->calculate_interactions(topo, conf, sim,
						 tid, m_omp_num_threads);
  }
  
#else
  
  // have to do all from here (probably it's only one, could unite this,
  // but then maybe it's clearer like it is...
  
  for(int i=0; i < m_omp_num_threads; ++i){
    m_nonbonded_set[i]->calculate_interactions(topo, conf, sim, 
					       i, m_omp_num_threads);
  }

#endif

  // add the forces, energies, virial...
  it = m_nonbonded_set.begin();

  const double ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & e = conf.current().energies;
  
  for( ; it != to; ++it){
    DEBUG(7, "adding forces from set " << it - m_nonbonded_set.begin());
    for(unsigned int i=0; i<topo.num_atoms(); ++i)
      conf.current().force(i) += (*it)->shortrange_storage().force(i);

    DEBUG(7, "adding energies from set " << it - m_nonbonded_set.begin());
    for(unsigned int i = 0; i < ljs; ++i){
      for(unsigned int j = 0; j < ljs; ++j){

	DEBUG(8, "set " << it - m_nonbonded_set.begin() << " group["
	      << i << "][" << j <<"] e_lj = " 
	      << (*it)->shortrange_storage().energies.lj_energy[i][j]
	      << " e_crf = " << (*it)->shortrange_storage().energies.crf_energy[i][j]);
      
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

    DEBUG(7, "\tadd perturbed energy derivatives");
    
    it = m_nonbonded_set.begin();
      
    for( ; it != to; ++it){
	
      configuration::Energy & pe = conf.current().perturbed_energy_derivatives;
	
      for(unsigned int i = 0; i < ljs; ++i){
	for(unsigned int j = 0; j < ljs; ++j){
      
	  assert(pe.lj_energy.size() > i);
	  assert(pe.lj_energy[i].size() > j);

	  assert((*it)->shortrange_storage().perturbed_energy_derivatives.lj_energy.size() > i);
	  assert((*it)->shortrange_storage().perturbed_energy_derivatives.lj_energy[i].size() > j);

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
 * calculate the interaction for a given atom pair.
 * SLOW! as it has to create the periodicity...
 */
int interaction::Nonbonded_Interaction::calculate_interaction
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 unsigned int atom_i, unsigned int atom_j,
 math::Vec & force, 
 double &e_lj, double &e_crf
 )
{
  assert(m_nonbonded_set.size() >= 1);
  return m_nonbonded_set[0]->calculate_interaction(topo, conf, sim,
						   atom_i, atom_j,
						   force, e_lj, e_crf);
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
int interaction::Nonbonded_Interaction::init(topology::Topology & topo,
					     configuration::Configuration & conf,
					     simulation::Simulation & sim,
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

    for(int i=0; i<m_omp_num_threads; ++i)
      m_nonbonded_set.push_back(new Perturbed_Nonbonded_Set(*m_pairlist_algorithm,
							    m_parameter));
    
  }
  else{
    for(int i=0; i<m_omp_num_threads; ++i)
      m_nonbonded_set.push_back(new Nonbonded_Set(*m_pairlist_algorithm, 
						  m_parameter));
  }
  
  std::vector<Nonbonded_Set *>::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();
  
  for( ; it != to; ++it){
    (*it)->init(topo, conf, sim, quiet);
  }

  check_spc_loop(topo, conf, sim, quiet);

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
     << std::setw(20) << "not measured: too expensive" << "\n";

  m_pairlist_algorithm->print_timing(os);

  os << "            "
     << std::setw(32) << std::left << "pairlist"
     << std::setw(20) 
     << m_pairlist_algorithm->timing() - m_longrange_timing<< "\n"
     << "\n";
}

void interaction::Nonbonded_Interaction::check_spc_loop
(
 topology::Topology const & topo,
 configuration::Configuration const & conf,
 simulation::Simulation & sim,
 bool quiet)
{
  DEBUG(7, "checking for spc interaction loops");
  
  if (sim.param().force.spc_loop == -1){

    sim.param().force.spc_loop = 0;
    if (!quiet)
      std::cout << "\tusing standard solvent loops (user request)\n";
    return;
  }
  
  if (topo.num_solvents() != 1 || topo.num_solvent_atoms(0) / topo.num_solvent_molecules(0) != 3 ||
      topo.num_solvent_molecules(0) < 1){

    sim.param().force.spc_loop = 0;
    if (!quiet)
	if (topo.num_solvents() > 0)
      std::cout << "\tusing standard solvent loops (num solvents doesn't match)\n"
		<< "\t\tnum solvents: " << topo.num_solvents() << "\n"
		<< "\t\tsolvent atoms: " << topo.num_solvent_atoms(0) / topo.num_solvent_molecules(0) << "\n"
		<< "\t\tmolecules: " << topo.num_solvent_molecules(0) << "\n\n";
	else
      std::cout << "\tusing standard solvent loops (no solvent present!)\n\n";
    return;
  }
  
  // check charges
  if (topo.charge()(topo.num_solute_atoms()) != -0.82 ||
      topo.charge()(topo.num_solute_atoms()+1) != 0.41 ||
      topo.charge()(topo.num_solute_atoms()+2) != 0.41){
    
    sim.param().force.spc_loop = 0;
    if (!quiet)
	std::cout << "\tusing standard solvent loops (charges don't match)\n"
		  << "\t\tO  : " << topo.charge()(topo.num_solute_atoms()) << "\n"
		  << "\t\tH1 : " << topo.charge()(topo.num_solute_atoms()+1) << "\n"
		  << "\t\tH2 : " << topo.charge()(topo.num_solute_atoms()+2) << "\n\n";
	  
    return;
  }
  
  // check lj parameters
  const lj_parameter_struct &lj_OO = 
    m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
			     topo.iac(topo.num_solute_atoms()));
  
  const lj_parameter_struct &lj_OH1 = 
    m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
			     topo.iac(topo.num_solute_atoms()+1));
  
  const lj_parameter_struct &lj_OH2 = 
    m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
			     topo.iac(topo.num_solute_atoms()+2));
  
  const lj_parameter_struct &lj_H1H2 = 
    m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()+1),
			     topo.iac(topo.num_solute_atoms()+2));
  
  if (lj_OO.c6 != 2.617346E-3 ||
      lj_OO.c12 != 2.634129E-6 ||
      lj_OH1.c6 != 0.0 ||
      lj_OH1.c12 != 0.0 ||
      lj_OH2.c6 != 0.0 ||
      lj_OH2.c12 != 0.0 ||
      lj_H1H2.c6 != 0.0 ||
      lj_H1H2.c12 != 0.0){
    
    sim.param().force.spc_loop = 0;
    if (!quiet)
      std::cout << "\tusing standard solvent loops (van der Waals parameter don't match)\n";
    return;
    
  }
  
  sim.param().force.spc_loop = 1;
  if (!quiet)
    std::cout << "\tusing spc solvent loops\n";
    return;
    
}
