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
::Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param)
  : m_pairlist_alg(pairlist_alg), 
    m_outerloop(param)
{
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::Nonbonded_Set
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 int tid, int num_threads)
{
  DEBUG(4, "Nonbonded_Set::calculate_interactions");
  
  // zero forces, energies, virial...
  m_shortrange_storage.zero();

  // need to update pairlist?
  if(!(sim.steps() % sim.param().pairlist.skip_step)){
    DEBUG(7, "\tdoing longrange...");
    
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
			  tid, topo.num_atoms(), num_threads);

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
  DEBUG(7, "\tshort range interactions");

  //  double shortrange_start = now();

  m_outerloop.lj_crf_outerloop(topo, conf, sim,
			       m_pairlist, m_shortrange_storage);
  
  // add 1,4 - interactions
  if (tid == 0){
    DEBUG(7, "\t1,4 - interactions");
    m_outerloop.one_four_outerloop(topo, conf, sim, m_shortrange_storage);
  
    // possibly do the RF contributions due to excluded atoms
    if(sim.param().longrange.rf_excluded){
      DEBUG(7, "\tRF excluded interactions and self term");
      m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_shortrange_storage);
    }
  }
  
  // add long-range force
  DEBUG(7, "\t(set) add long range forces");

  m_shortrange_storage.force += m_longrange_storage.force;
  
  // and long-range energies
  DEBUG(7, "\t(set) add long range energies");
  const unsigned int lj_e_size = unsigned(m_shortrange_storage.energies.lj_energy.size());
  
  for(unsigned int i = 0; i < lj_e_size; ++i){
    for(unsigned int j = 0; j < lj_e_size; ++j){
      m_shortrange_storage.energies.lj_energy[i][j] += 
	m_longrange_storage.energies.lj_energy[i][j];
      m_shortrange_storage.energies.crf_energy[i][j] += 
	m_longrange_storage.energies.crf_energy[i][j];
    }
  }

  // add longrange virial
  if (sim.param().pcouple.virial){
    DEBUG(7, "\t(set) add long range virial");
    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){

	DEBUG(8, "longrange virial = " << m_longrange_storage.virial_tensor(i,j)
	      << "\tshortrange virial = " << m_shortrange_storage.virial_tensor(i,j));

	m_shortrange_storage.virial_tensor(i,j) +=
	  m_longrange_storage.virial_tensor(i,j);
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
       bool quiet)
{
  // ?????
  // m_outerloop.initialize(sim);

  m_shortrange_storage.force.resize(conf.current().force.size());
  m_longrange_storage.force.resize(conf.current().force.size());

  m_shortrange_storage.energies.
    resize(unsigned(conf.current().energies.bond_energy.size()),
	   unsigned(conf.current().energies.kinetic_energy.size()));
  m_longrange_storage.energies.
    resize(unsigned(conf.current().energies.bond_energy.size()),
	   unsigned(conf.current().energies.kinetic_energy.size()));
  
  /*
  m_shortrange_storage.perturbed_energy_derivatives.resize
    (unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size()),
     unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()));

  m_longrange_storage.perturbed_energy_derivatives.resize
    (unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size()),
     unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()));
  */

  // and the pairlists
  pairlist().resize(topo.num_atoms());

  // check if we can guess the number of pairs
  const double vol = math::volume(conf.current().box, conf.boundary_type);
  if (vol){
    const double c3 = sim.param().pairlist.cutoff_short *
      sim.param().pairlist.cutoff_short *
      sim.param().pairlist.cutoff_short;
    
    const unsigned int pairs = 
      int(1.3 * topo.num_atoms() / vol * 4.0 / 3.0 * math::Pi * c3);

    if (!quiet)
      std::cout << "\n\testimated pairlist size (per atom) : "
		<< pairs << "\n\n";
    
    for(unsigned int i=0; i<topo.num_atoms(); ++i)
      pairlist()[i].reserve(pairs);
    
  }

  return 0;
}

