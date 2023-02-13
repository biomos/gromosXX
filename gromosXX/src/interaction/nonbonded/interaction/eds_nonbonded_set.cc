/**
 * @file eds_nonbonded_set.cc
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction_types.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"

#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_outerloop.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_term.h"

#include "../../../math/volume.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_set.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_set.h"

#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Eds_Nonbonded_Set
::Eds_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param, 
			  int rank, int num_threads)
  : Nonbonded_Set(pairlist_alg, param, rank, num_threads),
    m_eds_outerloop(param)   
{
    DEBUG(10, "EDS Nonbonded Set Constructor");
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::Eds_Nonbonded_Set
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "EDS Nonbonded_Set::calculate_interactions");
  
  // zero forces, energies, virial...
  m_storage.zero();

  // need to update pairlist?
  const bool pairlist_update = !(sim.steps() % sim.param().pairlist.skip_step);
  if(pairlist_update){
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
    
    if (topo.eds_perturbed_solute().atoms().size() > 0){
      m_pairlist_alg.update_perturbed(topo, conf, sim, 
				      pairlist(), perturbed_pairlist(),
				      m_rank, topo.num_atoms(), m_num_threads);
    }
    else{
      // assure it's empty
      assert(perturbed_pairlist().size() == topo.num_atoms());
      perturbed_pairlist().clear();
      
      
      m_pairlist_alg.update(topo, conf, sim, 
			    pairlist(),
			    m_rank, topo.num_atoms(), m_num_threads);
    }
    
  }  

  if(pairlist_update){
    //start_subtimer("longrange interactions (total)"); 

    m_outerloop.lj_crf_outerloop(topo, conf, sim,
			       m_pairlist.solute_long, m_pairlist.solvent_long,
                               m_longrange_storage, true /*longrange!*/, m_pairlist_alg.timer(), m_rank == 0);
     
    if (topo.eds_perturbed_solute().atoms().size() > 0){
      DEBUG(6, "\teds-perturbed long range");
      
      m_eds_outerloop.eds_lj_crf_outerloop(topo, conf, sim,
                                                 m_perturbed_pairlist.solute_long,
                                                 m_longrange_storage);
    }
    //stop_subtimer("longrange interactions (total)"); 
  }

  // calculate forces / energies
  DEBUG(6, "\tshort range interactions");
  //start_subtimer("shortrange interactions (total)");
  m_outerloop.lj_crf_outerloop(topo, conf, sim,
			       m_pairlist.solute_short, m_pairlist.solvent_short,
                               m_storage, false, m_pairlist_alg.timer(), m_rank == 0);
  
  if (topo.eds_perturbed_solute().atoms().size() > 0){
    DEBUG(6, "\tperturbed short range");
    m_eds_outerloop.eds_lj_crf_outerloop(topo, conf, sim, 
						     m_perturbed_pairlist.solute_short,
						     m_storage);
  }
  //stop_subtimer("shortrange interactions (total)");

  DEBUG(6, "\t1,4 - interactions");
  start_subtimer("1,4 interaction");
  m_outerloop.one_four_outerloop(topo, conf, sim, m_storage, m_rank, m_num_threads);
  stop_subtimer("1,4 interaction");

  start_subtimer("LJ exceptions");
  m_outerloop.lj_exception_outerloop(topo, conf, sim, m_storage,  m_rank, m_num_threads);
  stop_subtimer("LJ exceptions");

  if (sim.param().nonbonded.rf_excluded) {
    start_subtimer("RF excluded");
    DEBUG(6, "\tRF excluded interactions and self term");
    m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_storage, m_rank, m_num_threads);
    stop_subtimer("RF excluded");
  }

  // add 1,4 - interactions
  if (m_rank == 0) {
    if (topo.eds_perturbed_solute().atoms().size() > 0) {
      start_subtimer("1,4 interaction");
      DEBUG(6, "\teds-perturbed 1,4 - interactions");
      m_eds_outerloop.eds_one_four_outerloop(topo, conf, sim, m_storage);
      stop_subtimer("1,4 interaction");
    }
    

    // possibly do the RF contributions due to excluded atoms
    if (sim.param().nonbonded.rf_excluded) {
      start_subtimer("RF excluded");
      DEBUG(6, "\tRF excluded interactions and self term");
      if (topo.eds_perturbed_solute().atoms().size() > 0) {
        DEBUG(6, "\tperturbed RF excluded interactions and self term");
        m_eds_outerloop.eds_RF_excluded_outerloop(topo, conf, sim, m_storage);

      }
      stop_subtimer("RF excluded");
    }
  }
  
  // add longrange for EDS 
  const unsigned int numstates = m_storage.force_endstates.size();
  assert(m_storage.virial_tensor_endstates.size() == numstates &&
         m_longrange_storage.force_endstates.size() == numstates &&
         m_longrange_storage.virial_tensor_endstates.size() == numstates);
          
  // add long-range force
  DEBUG(6, "\t(set) add long range forces");
  m_storage.force += m_longrange_storage.force;
  // add eds long-range forces
  for(unsigned int i = 0; i < numstates; ++i){
    assert(m_storage.force_endstates[i].size()==m_longrange_storage.force_endstates[i].size());
    m_storage.force_endstates[i] += m_longrange_storage.force_endstates[i];
  }
   
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
  // and long-range energies of end states
  assert(m_storage.energies.eds_vi.size()==m_longrange_storage.energies.eds_vi.size());
  for(unsigned int i = 0; i < numstates; ++i){
    m_storage.energies.eds_vi[i] += m_longrange_storage.energies.eds_vi[i];
  }
  
  // add longrange virial
  if (sim.param().pcouple.virial){
    DEBUG(7, "\t(set) add long range virial");
    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){

	DEBUG(8, "longrange virial = " << m_longrange_storage.virial_tensor(i,j)
	      << "\tshortrange virial = " << m_storage.virial_tensor(i,j));

	m_storage.virial_tensor(i, j) +=
        m_longrange_storage.virial_tensor(i, j);
        
        // longrange virial of eds end states
        for(unsigned int k = 0; k < numstates; ++k){
          m_storage.virial_tensor_endstates[k](i, j) +=
          m_longrange_storage.virial_tensor_endstates[k](i, j);
        }
      }
    }
  }
  
  DEBUG(7, "(set) calculate interactions done!");

  return 0;
}

int interaction::Eds_Nonbonded_Set::update_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim)
{
  Nonbonded_Set::update_configuration(topo, conf, sim);
  
  // number of eds states
  const unsigned int numstates = m_storage.force_endstates.size();
  assert(m_storage.virial_tensor_endstates.size() == numstates);
  assert(conf.special().eds.force_endstates.size() == numstates);
  assert(conf.special().eds.virial_tensor_endstates.size() == numstates);
  assert(conf.current().energies.eds_vi.size() == numstates);
  assert(m_storage.energies.eds_vi.size() == numstates);

  // add storage forces to configuration
  for(unsigned int i = 0; i < numstates; ++i){
    for(unsigned int j = 0; j < topo.num_atoms(); ++j){
      conf.special().eds.force_endstates[i](j) += m_storage.force_endstates[i](j);
    }
  }
  
  // add storage energies to configuration
  for(unsigned int i = 0; i < numstates; ++i){
    conf.current().energies.eds_vi[i] += m_storage.energies.eds_vi[i];
    conf.current().energies.eds_vi_shift_extra_orig[i] += m_storage.energies.eds_vi_shift_extra_orig[i];
    conf.current().energies.eds_vi_shift_extra_phys[i] += m_storage.energies.eds_vi_shift_extra_phys[i];
  }
  
  // add storage virial 
  if (sim.param().pcouple.virial){
    DEBUG(7, "\tadd set virial");
    
    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){
        for(unsigned int k = 0; k < numstates; ++k){
          conf.special().eds.virial_tensor_endstates[k](i, j) +=
          m_storage.virial_tensor_endstates[k](i, j);
        }
      }
    }
  } // virial?

  return 0;
}

int interaction::Eds_Nonbonded_Set
::init(topology::Topology const & topo,
       configuration::Configuration const & conf,
       simulation::Simulation const & sim,
       std::ostream & os,
       bool quiet)
{
  DEBUG(7, "EDS Nonbonded Set :: init");

  Nonbonded_Set::init(topo, conf, sim, os, quiet); 
  perturbed_pairlist().resize(topo.num_atoms());
  
  m_storage.force_endstates.resize(sim.param().eds.numstates);
  m_storage.virial_tensor_endstates.resize(sim.param().eds.numstates);
  m_storage.energies.eds_vi.resize(sim.param().eds.numstates);
  m_storage.energies.eds_vi_shift_extra_orig.resize(sim.param().eds.numstates);
  m_storage.energies.eds_vi_shift_extra_phys.resize(sim.param().eds.numstates);

  m_longrange_storage.force_endstates.resize(sim.param().eds.numstates);
  m_longrange_storage.virial_tensor_endstates.resize(sim.param().eds.numstates);
  m_longrange_storage.energies.eds_vi.resize(sim.param().eds.numstates);
  
  assert(m_storage.force_endstates.size() 
      == m_longrange_storage.force_endstates.size());
  for(unsigned int i = 0; i < m_storage.force_endstates.size(); i++){
    m_storage.force_endstates[i].resize(topo.num_atoms());
    m_longrange_storage.force_endstates[i].resize(topo.num_atoms());
  }
  
  return 0;
}

