/**
 * @file perturbed_nonbonded_set.cc
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

#include <math/periodicity.h>

#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_term.h>

#include <interaction/nonbonded/interaction/perturbed_nonbonded_pair.h>

#include <math/volume.h>

#include <interaction/nonbonded/interaction/nonbonded_set.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_set.h>

#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Perturbed_Nonbonded_Set
::Perturbed_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param, 
			  int rank, int num_threads)
  : Nonbonded_Set(pairlist_alg, param, rank, num_threads),
    m_perturbed_outerloop(param),
    m_perturbed_pair(param)
{
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::Perturbed_Nonbonded_Set
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "Nonbonded_Set::calculate_interactions");

  const double l = topo.lambda();
  
  if (sim.param().perturbation.scaling){
    DEBUG(6, "scaling preparation");
    // calculate lambda primes and d lambda prime / d lambda derivatives
    std::map<std::pair<int, int>, std::pair<int, double> >::const_iterator
      it = topo.energy_group_lambdadep().begin(),
      to = topo.energy_group_lambdadep().end();
    
    for(unsigned int i=0; it != to; ++i, ++it){
      
      const double alpha = it->second.second;
      double lp = alpha * l * l + (1-alpha) * l;
      double dlp = (2 * l - 1.0) * alpha + 1;
      
      // some additional flexibility
      if (lp > 1.0) {
	lp = 1.0;
	dlp = 0.0;
      }
      
      else if (lp < 0.0){
	lp = 0.0;
	dlp = 0.0;
      }
      
      // -1 or not???
      assert(int(topo.lambda_prime().size()) > it->second.first);
      assert(int(topo.lambda_prime_derivative().size()) > it->second.first);
      assert(it->second.first >= 0);
      
      topo.lambda_prime()[it->second.first] = lp;
      topo.lambda_prime_derivative()[it->second.first] = dlp;
      
    }
  }
    
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
    DEBUG(6, "create a pairlist");
    
    if (topo.perturbed_solute().atoms().size() > 0){
      m_pairlist_alg.update_perturbed(topo, conf, sim, 
				      longrange_storage(),
				      pairlist(), perturbed_pairlist(),
				      m_rank, topo.num_atoms(), m_num_threads);
    }
    else{
      // assure it's empty
      assert(perturbed_pairlist().size() == topo.num_atoms());
      for(unsigned int i=0; i<topo.num_atoms(); ++i){
	perturbed_pairlist()[i].clear();
      }
      
      m_pairlist_alg.update(topo, conf, sim, 
			    longrange_storage(), pairlist(),
			    m_rank, topo.num_atoms(), m_num_threads);
    }
    
  }

  // calculate forces / energies
  DEBUG(6, "\tshort range interactions");

  m_outerloop.lj_crf_outerloop(topo, conf, sim,
			       m_pairlist, m_storage);

  if (topo.perturbed_solute().atoms().size() > 0){
    DEBUG(6, "\tperturbed short range");
    m_perturbed_outerloop.perturbed_lj_crf_outerloop(topo, conf, sim, 
						     m_perturbed_pairlist,
						     m_storage);
  }

  // add 1,4 - interactions
  if (m_rank == 0){
    DEBUG(6, "\t1,4 - interactions");
    m_outerloop.one_four_outerloop(topo, conf, sim, m_storage);

    if (topo.perturbed_solute().atoms().size() > 0){
      DEBUG(6, "\tperturbed 1,4 - interactions");
      m_perturbed_outerloop.perturbed_one_four_outerloop(topo, conf, sim, 
							 m_storage);
    }
    
    // possibly do the RF contributions due to excluded atoms
    if(sim.param().longrange.rf_excluded){
      DEBUG(6, "\tRF excluded interactions and self term");
      m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_storage);

      if (topo.perturbed_solute().atoms().size() > 0){
	DEBUG(6, "\tperturbed RF excluded interactions and self term");
	m_perturbed_outerloop.perturbed_RF_excluded_outerloop(topo, conf, sim,
							      m_storage);
      }
    }

    DEBUG(6, "\tperturbed pairs");
    m_perturbed_pair.perturbed_pair_outerloop(topo, conf, sim, m_storage);
    
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
    DEBUG(7, "\t(set) add long range virial");
    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){

	DEBUG(8, "longrange virial = " << m_longrange_storage.virial_tensor(i,j)
	      << "\tshortrange virial = " << m_storage.virial_tensor(i,j));

	m_storage.virial_tensor(i,j) +=
	  m_longrange_storage.virial_tensor(i,j);
      }
    }
  }
  
  // and long-range energy lambda-derivatives
  DEBUG(7, "(set) add long-range lambda-derivatives");
  
  const unsigned int lj_size 
    = unsigned(m_storage.perturbed_energy_derivatives.lj_energy.size());
  
  for(unsigned int i = 0; i < lj_size; ++i){
    for(unsigned int j = 0; j < lj_size; ++j){
      
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy.size() > i);
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy[i].size() > j);
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy.size() > j);
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy[j].size() > i);
      assert(m_longrange_storage.perturbed_energy_derivatives.
	     lj_energy.size() > i);
      assert(m_longrange_storage.perturbed_energy_derivatives.
	     lj_energy[i].size() > j);
      assert(m_longrange_storage.perturbed_energy_derivatives.
	     lj_energy.size() > j);
      assert(m_longrange_storage.perturbed_energy_derivatives.
	     lj_energy[j].size() > i);
      
      m_storage.perturbed_energy_derivatives.lj_energy[i][j] += 
	m_longrange_storage.perturbed_energy_derivatives.lj_energy[i][j];
      
      m_storage.perturbed_energy_derivatives.crf_energy[i][j] += 
	m_longrange_storage.perturbed_energy_derivatives.crf_energy[i][j];
      
    }
  }
  
  DEBUG(7, "(set) calculate interactions done!");

  return 0;
}

int interaction::Perturbed_Nonbonded_Set::update_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim)
{
  const int ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & pe = conf.current().perturbed_energy_derivatives;

  Nonbonded_Set::update_configuration(topo, conf, sim);
	
  for(int i = 0; i < ljs; ++i){
    for(int j = 0; j < ljs; ++j){
      
      assert(pe.lj_energy.size() > unsigned(i));
      assert(pe.lj_energy[i].size() > unsigned(j));

      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy.size() > unsigned(i));
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy[i].size() > (unsigned(j)));
	  
      pe.lj_energy[i][j] += 
	m_storage.perturbed_energy_derivatives.
	lj_energy[i][j];
      pe.crf_energy[i][j] += 
	m_storage.perturbed_energy_derivatives.
	crf_energy[i][j];
    }
  }
  return 0;
}

/**
 * calculate the hessian for a given atom.
 * this will be VERY SLOW !
 */
int interaction::Perturbed_Nonbonded_Set
::calculate_hessian(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    unsigned int atom_i, unsigned int atom_j,
		    math::Matrix & hessian){
  
  if (topo.is_perturbed(atom_i) ||
      topo.is_perturbed(atom_j)){
    assert(false);
    return -1;
  }

  return Nonbonded_Set::calculate_hessian(topo, conf, sim,
					  atom_i, atom_j, hessian);
}

int interaction::Perturbed_Nonbonded_Set
::init(topology::Topology const & topo,
       configuration::Configuration const & conf,
       simulation::Simulation const & sim,
       std::ostream & os,
       bool quiet)
{
  DEBUG(7, "Perturbed Nonbonded Set :: init");
  DEBUG(7, "sizes: " 
	<< unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size())
	<< " - " 
	<< unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()));

  Nonbonded_Set::init(topo, conf, sim, os, quiet);
  
  m_storage.perturbed_energy_derivatives.resize
    (unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size()),
     unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()));

  m_longrange_storage.perturbed_energy_derivatives.resize
    (unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size()),
     unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()));

  perturbed_pairlist().resize(topo.num_atoms());

  return 0;
}

