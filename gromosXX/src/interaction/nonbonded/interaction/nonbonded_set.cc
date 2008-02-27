/**
 * @file nonbonded_set.cc
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>

#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>

#include <math/periodicity.h>
#include <math/volume.h>

#include <interaction/nonbonded/interaction/nonbonded_set.h>

#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Nonbonded_Set
::Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param,
		int rank, int num_threads)
  : Nonbonded_Set_Interface(pairlist_alg, param, rank, num_threads),
    m_outerloop(param)
{
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::Nonbonded_Set
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "Nonbonded_Set::calculate_interactions");
  
  // zero forces, energies, virial...
  m_storage.zero();

  // need to update pairlist?
  if(!(sim.steps() % sim.param().pairlist.skip_step)){
    DEBUG(6, "\tdoing longrange...");
    
    //====================
    // create a pairlist
    //====================
    
    // zero the longrange forces, energies, virial
    m_longrange_storage.zero();

    // parallelisation using STRIDE:
    // chargegroup based pairlist can only use this one!!!!
    // TODO:
    // move decision to pairlist!!!
    m_pairlist_alg.update(topo, conf, sim, 
			  longrange_storage(), pairlist(),
			  m_rank, topo.num_atoms(), m_num_threads);

    /*
    sleep(2*tid);
    
    std::cout << "PRINTING OUT THE PAIRLIST\n\n";
    for(unsigned int i=0; i<100; ++i){
      if (i >= pairlist().size()) break;

      std::cout << "\n\n--------------------------------------------------";
      std::cout << "\n" << i;
      for(unsigned int j=0; j<pairlist()[i].size(); ++j){

	if (j % 10 == 0) std::cout << "\n\t";
	std::cout << std::setw(7) << pairlist()[i][j];
      }
    }
    */
  }

  // calculate forces / energies
  DEBUG(6, "\tshort range interactions");
  
  m_outerloop.lj_crf_outerloop(topo, conf, sim,
			       m_pairlist, m_storage);
  
  // add 1,4 - interactions
  if (m_rank == 0){
    DEBUG(6, "\t1,4 - interactions");
    m_outerloop.one_four_outerloop(topo, conf, sim, m_storage);
  
    // possibly do the RF contributions due to excluded atoms
    if(sim.param().longrange.rf_excluded){
      DEBUG(7, "\tRF excluded interactions and self term");
      m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_storage);
    }
  }
  
  // add long-range force
  DEBUG(6, "\t(set) add long range forces");

  m_storage.force += m_longrange_storage.force;
  
  // and long-range energies
  DEBUG(6, "\t(set) add long range energies");
  const unsigned int lj_e_size = unsigned(m_storage.energies.lj_energy.size());
  
  for(unsigned int i = 0; i < lj_e_size; ++i){
    for(unsigned int j = 0; j < lj_e_size; ++j){
      m_storage.energies.lj_energy[i][j] += 
	m_longrange_storage.energies.lj_energy[i][j];
      m_storage.energies.crf_energy[i][j] += 
	m_longrange_storage.energies.crf_energy[i][j];
    }
  }

  // add longrange virial
  if (sim.param().pcouple.virial){
    DEBUG(6, "\t(set) add long range virial");
    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){

	DEBUG(8, "longrange virial = " << m_longrange_storage.virial_tensor(i,j)
	      << "\tshortrange virial = " << m_storage.virial_tensor(i,j));

	m_storage.virial_tensor(i,j) +=
	  m_longrange_storage.virial_tensor(i,j);
      }
    }
  }
  
  return 0;
}

int interaction::Nonbonded_Set::update_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim)
{
  const int ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & e = conf.current().energies;

  // use the IMPULSE method for multiple time stepping
  if (sim.param().multistep.steps > 1){
    int steps = sim.param().multistep.steps;
    if (sim.param().multistep.boost == 0)
      steps = 1;
    
    // only add when calculated
    if ((sim.steps() % steps) == 0){

      // std::cerr << "\tadding boosted (" << steps << ") non-bonded forces" << std::endl;

      for(unsigned int i=0; i<topo.num_atoms(); ++i)
	conf.current().force(i) += steps * m_storage.force(i);
    }
    else{
      // std::cerr << "\tnot adding non-bonded forces" << std::endl;
    }
    
  }
  else{
    for(unsigned int i=0; i<topo.num_atoms(); ++i)
      conf.current().force(i) += m_storage.force(i);
  }
  
  // (MULTISTEP: and keep energy constant)
  for(int i = 0; i < ljs; ++i){
    for(int j = 0; j < ljs; ++j){
      
      e.lj_energy[i][j] += 
	m_storage.energies.lj_energy[i][j];
      e.crf_energy[i][j] += 
	m_storage.energies.crf_energy[i][j];
    }
  }

  // (MULTISTEP: and the virial???)
  if (sim.param().pcouple.virial){
    DEBUG(7, "\tadd set virial");

    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){

	conf.current().virial_tensor(i,j) +=
	  m_storage.virial_tensor(i,j);
      }
    }
  }
  return 0;
}

/**
 * calculate the interaction for a given atom pair.
 * SLOW! as it has to create the periodicity...
 */
int interaction::Nonbonded_Set::calculate_interaction
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 unsigned int atom_i, unsigned int atom_j,
 math::Vec & force, 
 double &e_lj, double &e_crf
 )
{
  return m_outerloop.calculate_interaction(topo, conf, sim,
					   atom_i, atom_j,
					   force, e_lj, e_crf);
}


/**
 * calculate the hessian for a given atom.
 * this will be VERY SLOW !
 */
int interaction::Nonbonded_Set
::calculate_hessian(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    unsigned int atom_i, unsigned int atom_j,
		    math::Matrix & hessian){
  
  return m_outerloop.calculate_hessian(topo, conf, sim,
				       atom_i, atom_j, hessian,
				       m_pairlist);
}

int interaction::Nonbonded_Set
::init(topology::Topology const & topo,
       configuration::Configuration const & conf,
       simulation::Simulation const & sim,
       std::ostream & os,
       bool quiet)
{
  // ?????
  // m_outerloop.initialize(sim);

  // std::cerr << "nonbonded set: init" << std::endl;
  
  const int num_atoms = topo.num_atoms();

  m_storage.force.resize(num_atoms);
  m_longrange_storage.force.resize(num_atoms);

  m_storage.energies.
    resize(unsigned(conf.current().energies.bond_energy.size()),
	   unsigned(conf.current().energies.kinetic_energy.size()));
  m_longrange_storage.energies.
    resize(unsigned(conf.current().energies.bond_energy.size()),
	   unsigned(conf.current().energies.kinetic_energy.size()));
  
  // and the pairlists
  DEBUG(10, "pairlist size: " << num_atoms);
  pairlist().resize(num_atoms);

  // check if we can guess the number of pairs
  const double vol = math::volume(conf.current().box, conf.boundary_type);
  if (vol){
    const double c3 = sim.param().pairlist.cutoff_short *
      sim.param().pairlist.cutoff_short *
      sim.param().pairlist.cutoff_short;
    
    const unsigned int pairs = 
      int(1.3 * num_atoms / vol * 4.0 / 3.0 * math::Pi * c3);

    if (!quiet)
      os << "\n\testimated pairlist size (per atom) : "
	 << pairs << "\n";
    
    for(int i=0; i<num_atoms; ++i)
      pairlist()[i].reserve(pairs);
  }
  return 0;
}

