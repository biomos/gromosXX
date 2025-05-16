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

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_pair.h"

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
    m_eds_outerloop(param), m_perturbed_pair(param)  
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
    
    DEBUG(10, "\teds_perturbed_solute().atoms().size() " << topo.eds_perturbed_solute().atoms().size());
    DEBUG(10, "\tperturbed_solute().atoms().size() " << topo.perturbed_solute().atoms().size());
    if (topo.eds_perturbed_solute().atoms().size() +topo.perturbed_solute().atoms().size() > 0){
      DEBUG(10, "update perturbed");
      
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

  if (sim.param().polarise.cos) {
    //===============
    // polarisation
    //===============
    // Only compatible with TI for now!

    // calculate explicit polarisation of the molecules
    DEBUG(6, "\texplicit polarisation");
    start_timer("explicit polarisation");
    if (topo.perturbed_solute().atoms().size() > 0) {
      m_eds_outerloop.perturbed_electric_field_outerloop(topo, conf, sim,
                                       m_pairlist, m_perturbed_pairlist,
				       m_storage, m_longrange_storage, m_rank);
    } else {
      m_outerloop.electric_field_outerloop(topo, conf, sim, m_pairlist, 
				       m_storage, m_longrange_storage, m_rank);
    }
    stop_timer("explicit polarisation");
  }   

  if(pairlist_update){
    start_timer("longrange");
    DEBUG(6, "\tlong range interactions");
    
    m_outerloop.lj_crf_outerloop(topo, conf, sim,
			       m_pairlist.solute_long, m_pairlist.solvent_long,
                               m_longrange_storage, true /*longrange!*/, m_pairlist_alg.timer(), m_rank == 0);

    DEBUG(9,"nonEDS: virial longrange: " );
    for (unsigned int a=0; a<3; ++a){ 
      for (unsigned int b=0; b<3; ++b){
        DEBUG(9,"\t"<< m_longrange_storage.virial_tensor(a,b));
      }
    }
    
    if (topo.eds_perturbed_solute().atoms().size() + topo.perturbed_solute().atoms().size() > 0){
      DEBUG(6, "\teds-perturbed long range");
      
      m_eds_outerloop.eds_lj_crf_outerloop(topo, conf, sim,
                                                 m_perturbed_pairlist.solute_long,
                                                 m_longrange_storage);
    }
    stop_timer("longrange");
  }

  // calculate forces / energies
   DEBUG(6, "\tshort range interactions");
  start_timer("shortrange");
  m_outerloop.lj_crf_outerloop(topo, conf, sim,
			       m_pairlist.solute_short, m_pairlist.solvent_short,
                               m_storage, false, m_pairlist_alg.timer(), m_rank == 0);
  DEBUG(9,"nonEDS: virial shortrange: " );
  for (unsigned int a=0; a<3; ++a){ 
    for (unsigned int b=0; b<3; ++b){
      DEBUG(9,"\t"<< m_storage.virial_tensor(a,b));
    }
  }
  
  if (topo.eds_perturbed_solute().atoms().size() + topo.perturbed_solute().atoms().size() > 0){
    DEBUG(6, "\teds-perturbed short range");
    m_eds_outerloop.eds_lj_crf_outerloop(topo, conf, sim, 
					 m_perturbed_pairlist.solute_short,
					 m_storage);
  }
  stop_timer("shortrange");

  DEBUG(6, "\t1,4 - interactions");
  start_timer("1,4 interaction");
  m_outerloop.one_four_outerloop(topo, conf, sim, m_storage, m_rank, m_num_threads);
  stop_timer("1,4 interaction");

  start_timer("LJ exceptions");
  m_outerloop.lj_exception_outerloop(topo, conf, sim, m_storage,  m_rank, m_num_threads);
  stop_timer("LJ exceptions");

  if (sim.param().nonbonded.rf_excluded) {
    start_timer("RF excluded");
    DEBUG(6, "\tRF excluded interactions and self term");
    m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_storage, m_rank, m_num_threads);
    stop_timer("RF excluded");
  }

  // add 1,4 - interactions
  if (m_rank == 0){
    if (sim.param().polarise.cos) {
      // Polarisation is only commpatible with TI for now!!
      start_timer("polarisation self-energy");
      if (topo.perturbed_solute().atoms().size() > 0) {
        m_eds_outerloop.perturbed_self_energy_outerloop(topo, conf, sim, m_storage);
      } else {
        m_outerloop.self_energy_outerloop(topo, conf, sim, m_storage);
      }
      stop_timer("polarisation self-energy");
    }

    if (topo.eds_perturbed_solute().atoms().size() + topo.perturbed_solute().atoms().size() > 0) {
      start_timer("1,4 interaction");
      DEBUG(6, "\teds-perturbed 1,4 - interactions");
      m_eds_outerloop.eds_one_four_outerloop(topo, conf, sim, m_storage);
      stop_timer("1,4 interaction");
    }
    

    // possibly do the RF contributions due to excluded atoms
    if (sim.param().nonbonded.rf_excluded) {
      start_timer("RF excluded");
      DEBUG(6, "\tRF excluded interactions and self term");
      if (topo.eds_perturbed_solute().atoms().size() +  topo.perturbed_solute().atoms().size() > 0) {
        DEBUG(6, "\teds-perturbed RF excluded interactions and self term");
        m_eds_outerloop.eds_RF_excluded_outerloop(topo, conf, sim, m_storage);
      }
      stop_timer("RF excluded");
    }

    DEBUG(6, "\tperturbed pairs");
    start_timer("perturbed pairs");
    m_perturbed_pair.perturbed_pair_outerloop(topo, conf, sim, m_storage);
    stop_timer("perturbed pairs");
  }
  
  // add longrange for multiAEDS
  if (sim.param().eds.numsites > 1) {
    const unsigned int numentries = m_storage.force_mult_endstates.size();
    assert( m_storage.virial_tensor_mult_endstates.size() == numentries &&
            m_longrange_storage.force_mult_endstates.size() == numentries &&
            m_longrange_storage.virial_tensor_mult_endstates.size() == numentries);

    // add multiAEDS long-range forces
    for (auto i: m_longrange_storage.force_mult_endstates){
      m_storage.force_mult_endstates[i.first] += i.second;
      //assert(m_storage.force_endstates[i].size()==m_longrange_storage.force_endstates[i].size());
      //m_storage.force_endstates[i] += m_longrange_storage.force_endstates[i];
    }

    // add multiAEDS long-range energies
    for (auto i: m_longrange_storage.energies.eds_mult_vi){
      m_storage.energies.eds_mult_vi[i.first] += i.second;
    }

    // add multiAEDS longrange virial 
    if (sim.param().pcouple.virial){
      DEBUG(7, "\tmultiAEDS, add longrange virial contributions to shortrange in storage");
      for (unsigned int i=0; i<3; ++i){ //do we need this? or can we do += on the matrix?
        for (unsigned int j=0; j<3; ++j){
          for (auto k: m_longrange_storage.virial_tensor_mult_endstates){
            DEBUG(8, "longrange virial multiAEDS= " << k.second(i,j)
  	          << "\tshortrange virial = " << m_storage.virial_tensor_mult_endstates[k.first](i,j) 
              << "\ttotal: " << m_storage.virial_tensor_mult_endstates[k.first](i,j) + k.second(i,j) );
            m_storage.virial_tensor_mult_endstates[k.first](i, j) += k.second(i,j);
          }
        }
      }
    } //end virial


  } else {
    // add longrange for EDS 
    const unsigned int numstates = m_storage.force_endstates.size();
    assert(m_storage.virial_tensor_endstates.size() == numstates &&
          m_longrange_storage.force_endstates.size() == numstates &&
          m_longrange_storage.virial_tensor_endstates.size() == numstates);
            
    // add eds long-range forces
    for(unsigned int i = 0; i < numstates; ++i){
      assert(m_storage.force_endstates[i].size()==m_longrange_storage.force_endstates[i].size());
      m_storage.force_endstates[i] += m_longrange_storage.force_endstates[i];
    }
    
    // and long-range energies of end states
    assert(m_storage.energies.eds_vi.size()
        ==m_longrange_storage.energies.eds_vi.size());
    assert(m_storage.perturbed_energy_derivatives.eds_vi.size()
        ==m_longrange_storage.perturbed_energy_derivatives.eds_vi.size());
    assert(m_storage.perturbed_energy_derivatives.eds_vi.size() == numstates);
    DEBUG(6, "\t(set) add long-range energies of end states");
    for(unsigned int i = 0; i < numstates; ++i){
      m_storage.energies.eds_vi[i] += m_longrange_storage.energies.eds_vi[i];
      DEBUG(6, "\t\t(set) eds_vi " << m_longrange_storage.energies.eds_vi[i]);
      m_storage.perturbed_energy_derivatives.eds_vi[i] += m_longrange_storage.perturbed_energy_derivatives.eds_vi[i];
      DEBUG(6, "\t\t(set) eds_dvi " << m_longrange_storage.perturbed_energy_derivatives.eds_vi[i]);
    }
    
    // add longrange virial eds
    if (sim.param().pcouple.virial){
      DEBUG(7, "\t(set) add long range virial eds");
      for(unsigned int i=0; i<3; ++i){
        for(unsigned int j=0; j<3; ++j){
          for(unsigned int k = 0; k < numstates; ++k){
            m_storage.virial_tensor_endstates[k](i, j) +=
            m_longrange_storage.virial_tensor_endstates[k](i, j);
          }
        }
      }
    } //end virial
  }

    // add long-range force
    DEBUG(6, "\t(set) add long range forces");
    m_storage.force += m_longrange_storage.force;

  // and non-eds long-range energies
  DEBUG(6, "\t(set) add long range energies");
  const unsigned int lj_e_size = unsigned(m_storage.energies.lj_energy.size());
  
  for(unsigned int i = 0; i < lj_e_size; ++i){
    for(unsigned int j = 0; j < lj_e_size; ++j){
      m_storage.energies.lj_energy[i][j] += 
	m_longrange_storage.energies.lj_energy[i][j];
      m_storage.energies.crf_energy[i][j] += 
	m_longrange_storage.energies.crf_energy[i][j];
      m_storage.perturbed_energy_derivatives.lj_energy[i][j] +=
	m_longrange_storage.perturbed_energy_derivatives.lj_energy[i][j];
      m_storage.perturbed_energy_derivatives.crf_energy[i][j] +=
	m_longrange_storage.perturbed_energy_derivatives.crf_energy[i][j];
    }
  }

  // add longrange virial, non-eds
  if (sim.param().pcouple.virial){
    DEBUG(7, "\t(set) add long range virial, non-eds (eds_nonbonded_set)");
    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){

	      DEBUG(8, "longrange virial = " << m_longrange_storage.virial_tensor(i,j)
	         << "\tshortrange virial = " << m_storage.virial_tensor(i,j));

	      m_storage.virial_tensor(i, j) +=
	        m_longrange_storage.virial_tensor(i, j);
      }
    }
  }

  // Extended TI. Not compatible with EDS yet
  if (sim.param().precalclam.nr_lambdas){
    for(unsigned int i = 0; i < sim.param().precalclam.nr_lambdas; ++i){
      for(unsigned int j = 0; j < lj_e_size; ++j){
        for(unsigned int k = 0; k < lj_e_size; ++k){
          m_storage.energies.A_lj_energy[i][j][k] +=
            m_longrange_storage.energies.A_lj_energy[i][j][k];
          m_storage.energies.B_lj_energy[i][j][k] +=
            m_longrange_storage.energies.B_lj_energy[i][j][k];
          m_storage.energies.A_crf_energy[i][j][k] +=
            m_longrange_storage.energies.A_crf_energy[i][j][k];
          m_storage.energies.B_crf_energy[i][j][k] +=
            m_longrange_storage.energies.B_crf_energy[i][j][k];
          m_storage.perturbed_energy_derivatives.A_lj_energy[i][j][k] +=
            m_longrange_storage.perturbed_energy_derivatives.A_lj_energy[i][j][k];
          m_storage.perturbed_energy_derivatives.B_lj_energy[i][j][k] +=
            m_longrange_storage.perturbed_energy_derivatives.B_lj_energy[i][j][k];
          m_storage.perturbed_energy_derivatives.A_crf_energy[i][j][k] +=
            m_longrange_storage.perturbed_energy_derivatives.A_crf_energy[i][j][k];
          m_storage.perturbed_energy_derivatives.B_crf_energy[i][j][k] +=
            m_longrange_storage.perturbed_energy_derivatives.B_crf_energy[i][j][k];
        }
      }
    }
  } // Extended TI

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
  // this should, theoretically, also have taken care of all the stuff that the 
  // eds_outer_loop has written to m_storage 

  // multiAEDS
  if (sim.param().eds.numsites > 1) {
    const unsigned int numentries = m_storage.force_mult_endstates.size();
    assert(m_storage.virial_tensor_mult_endstates.size() == numentries);
    assert(conf.special().eds.force_mult_endstates.size() == numentries);
    assert(conf.special().eds.virial_tensor_mult_endstates.size() == numentries);
    assert(m_storage.energies.eds_mult_vi.size() == numentries);
    assert(conf.current().energies.eds_mult_vi.size() == numentries);

    double start_time = util::now();
    //add storage forces to configuration (separately for each atom?)
    //for (auto i: m_storage.force_mult_endstates){
    //  conf.special().eds.force_mult_endstates[i.first] += i.second;
    //}
    auto& force_mult_endstates = conf.special().eds.force_mult_endstates;
    DEBUG(7,"Number of elements force_mult_endstates: " << force_mult_endstates.size())
    for (const auto& [key, value] : m_storage.force_mult_endstates) {
      force_mult_endstates[key] += value;
    }
    DEBUG(7," Time for eds forces: " << util::now()-start_time );

    start_time = util::now();
    //add storage energies to configuration
    for (auto i: m_storage.energies.eds_mult_vi){
      conf.current().energies.eds_mult_vi[i.first] += i.second;
    }
    DEBUG(7," Time for eds energies: " << util::now()-start_time);

    start_time = util::now();
    // add storage virial 
    if (sim.param().pcouple.virial){
      DEBUG(7, "\tadd set virial for multiAEDS");
      DEBUG(9,"\tadding virial_tensor_mult_endstates from m_storage to conf.special().eds");
      for(unsigned int i=0; i<3; ++i){
        for(unsigned int j=0; j<3; ++j){
          for(auto k: m_storage.virial_tensor_mult_endstates){
            DEBUG(9,"\tk: " << k.first[0] << k.first[1] << k.first[2] << k.first[3] << " conf.special().eds.virial_tensor_mult_endstates= " 
            <<conf.special().eds.virial_tensor_mult_endstates[k.first](i,j) << "\tadding= " << k.second(i,j));
            conf.special().eds.virial_tensor_mult_endstates[k.first](i, j) +=
                k.second(i, j);
          }
        }
      }
    } 
    DEBUG(7," Time for eds virial: " << util::now()-start_time);

  } else {
    // number of eds states
    const unsigned int numstates = m_storage.force_endstates.size();
    assert(m_storage.virial_tensor_endstates.size() == numstates);
    assert(conf.special().eds.force_endstates.size() == numstates);
    assert(conf.special().eds.virial_tensor_endstates.size() == numstates);
    assert(conf.current().energies.eds_vi.size() == numstates);
    assert(conf.current().perturbed_energy_derivatives.eds_vi.size() == numstates);
    assert(m_storage.energies.eds_vi.size() == numstates);
    assert(m_storage.perturbed_energy_derivatives.eds_vi.size() == numstates);
    
    // add storage forces to configuration
    for(unsigned int i = 0; i < numstates; ++i){
      for(unsigned int j = 0; j < topo.num_atoms(); ++j){
        conf.special().eds.force_endstates[i](j) += m_storage.force_endstates[i](j);
      }
    }
    
    // add storage energies to configuration
    for(unsigned int i = 0; i < numstates; ++i){
      conf.current().energies.eds_vi[i] += m_storage.energies.eds_vi[i];
      conf.current().perturbed_energy_derivatives.eds_vi[i] += m_storage.perturbed_energy_derivatives.eds_vi[i]; 
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
    } 
  } //end eds 

  // we have written non-eds interactions directly to the forces, energies and derivates
  // The forces and energies have been taken care of by the call 
  // to Nonbonded_Set::update_configuration, but we still need to include the derivatives 
  const int ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & pe = conf.current().perturbed_energy_derivatives;
  
  // energy derivatives
  for (int i = 0; i < ljs; ++i) {
    for (int j = 0; j < ljs; ++j) {
      pe.lj_energy[i][j] += 
	        m_storage.perturbed_energy_derivatives.lj_energy[i][j];
      pe.crf_energy[i][j] += 
	        m_storage.perturbed_energy_derivatives.crf_energy[i][j];
    }
    pe.self_energy[i] += m_storage.perturbed_energy_derivatives.self_energy[i];
  }

  // Save extended TI derivatives
  const unsigned int nr_lambdas = unsigned(m_storage.energies.A_lj_total.size());

  for(unsigned int i=0; i < nr_lambdas; ++i) {
    for(int j=0; j < ljs; ++j) {
      for(int k=0; k < ljs; ++k) {
        pe.A_lj_energy[i][j][k] +=
            m_storage.perturbed_energy_derivatives.A_lj_energy[i][j][k];
        pe.B_lj_energy[i][j][k] +=
            m_storage.perturbed_energy_derivatives.B_lj_energy[i][j][k];
        pe.A_crf_energy[i][j][k] +=
            m_storage.perturbed_energy_derivatives.A_crf_energy[i][j][k];
        pe.B_crf_energy[i][j][k] +=
            m_storage.perturbed_energy_derivatives.B_crf_energy[i][j][k];
      }
    }
  } 
  
  return 0;
}


/**
 * calculate the hessian for a given atom.
 * this will be VERY SLOW !
 */
int interaction::Eds_Nonbonded_Set
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

  m_storage.perturbed_energy_derivatives.resize
    (unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size()),
     unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()),
     unsigned(sim.param().precalclam.nr_lambdas));

  m_longrange_storage.perturbed_energy_derivatives.resize
    (unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size()),
     unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()),
     unsigned(sim.param().precalclam.nr_lambdas));

  m_storage.force_endstates.resize(sim.param().eds.numstates);
  m_storage.virial_tensor_endstates.resize(sim.param().eds.numstates);
  m_storage.energies.eds_vi.resize(sim.param().eds.numstates);
  m_storage.perturbed_energy_derivatives.eds_vi.resize(sim.param().eds.numstates);
  
  m_longrange_storage.force_endstates.resize(sim.param().eds.numstates);
  m_longrange_storage.virial_tensor_endstates.resize(sim.param().eds.numstates);
  m_longrange_storage.energies.eds_vi.resize(sim.param().eds.numstates);
  m_longrange_storage.perturbed_energy_derivatives.eds_vi.resize(sim.param().eds.numstates);

  assert(m_storage.force_endstates.size() 
      == m_longrange_storage.force_endstates.size());
  for(unsigned int i = 0; i < m_storage.force_endstates.size(); i++){
    m_storage.force_endstates[i].resize(topo.num_atoms());
    m_longrange_storage.force_endstates[i].resize(topo.num_atoms());
  }

  //MULTIAEDS
  //initialize the maps of eds_mult_vi 
  //do we need this?
  //maybe only need to zero?
  for (auto i: sim.param().eds.site_state_pairs){
    m_storage.energies.eds_mult_vi[i];
    m_longrange_storage.energies.eds_mult_vi[i];
    m_storage.force_mult_endstates[i].resize(topo.num_atoms());
    m_longrange_storage.force_mult_endstates[i].resize(topo.num_atoms());
    m_storage.virial_tensor_mult_endstates[i];
    m_longrange_storage.virial_tensor_mult_endstates[i];
  }
  /*
  for( int site_i = 0; site_i < sim.param().eds.numsites; site_i++){
    for( int site_j = site_i; site_j < sim.param().eds.numsites; site_j++){
      for( int state_i =0; state_i < sim.param().eds.multnumstates[site_i]+1; state_i++){
        if( site_i == site_j){
          std::vector<int> states_i_j = {site_i, state_i, site_i, state_i};
	        m_storage.energies.eds_mult_vi[states_i_j];
          m_longrange_storage.energies.eds_mult_vi[states_i_j];
          m_storage.force_mult_endstates[states_i_j].resize(topo.num_atoms());//?
          m_longrange_storage.force_mult_endstates[states_i_j].resize(topo.num_atoms());//?
          m_storage.virial_tensor_mult_endstates[states_i_j];
          m_longrange_storage.virial_tensor_mult_endstates[states_i_j];
        }
        else{
          for( int state_j=0; state_j < sim.param().eds.multnumstates[site_j]; state_j++){
            std::vector<int> states_i_j = {site_i, state_i, site_j, state_j};
            m_storage.energies.eds_mult_vi[states_i_j];
            m_longrange_storage.energies.eds_mult_vi[states_i_j];
            m_storage.force_mult_endstates[states_i_j].resize(topo.num_atoms());//?
            m_longrange_storage.force_mult_endstates[states_i_j].resize(topo.num_atoms());//?
            m_storage.virial_tensor_mult_endstates[states_i_j];
            m_longrange_storage.virial_tensor_mult_endstates[states_i_j];
          }
        }
      }
    }
  }
  */
  return 0;
}

