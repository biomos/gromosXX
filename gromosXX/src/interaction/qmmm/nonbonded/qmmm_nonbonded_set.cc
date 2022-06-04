/**
 * @file qmmm_nonbonded_set.cc
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../simulation/parameter.h"

#include "../../../interaction/interaction.h"

#include "../../../util/debug.h"
#include "../../../util/timing.h"

#include "../../../configuration/energy.h"

#include "../../../math/volume.h"

#include "qm_zone.h"
#include "qmmm_nonbonded_set.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

/**
 * Constructor.
 */
interaction::QMMM_Nonbonded_Set::QMMM_Nonbonded_Set(QM_Zone& qm_zone
                                                  , util::Algorithm_Timer& timer
                                                  , Nonbonded_Parameter& param
		                                              , int rank, int num_threads)
  : m_qm_zone(qm_zone)
  , m_timer(timer)
  , m_outerloop(param)
  , m_rank(rank)
  , m_num_threads(num_threads)
{}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::QMMM_Nonbonded_Set::calculate_interactions(
                                        topology::Topology & topo
                                      , configuration::Configuration & conf
			                                , simulation::Simulation & sim) {
  DEBUG(4, "QMMM_Nonbonded_Set::calculate_interactions");
  
  // zero forces, energies, virial...
  m_storage.zero();

  // update the pairlist
  start_timer("pairlist update");
  m_qm_zone.update_pairlist(topo, sim, m_pairlist, m_rank, topo.num_atoms(), m_num_threads);
  stop_timer("pairlist update");

  // calculate forces / energies
  DEBUG(6, "\tQMMM LJ interactions");
  start_timer("nonbonded LJ");

  // With QMMM we do only single range cutoff - shortrange
  m_outerloop.lj_outerloop(topo, conf, sim,
          m_pairlist.solute_short, m_pairlist.solvent_short,
          m_storage, false, m_timer,
          m_rank == 0);

  stop_timer("nonbonded LJ");
  //Possibly do one_four_interaction and LJ_exception here

  DEBUG(6, "\tQMMM 1,4 - interactions");
  start_timer("1,4 interaction");
  m_outerloop.one_four_outerloop(topo, conf, sim, m_storage, m_rank, m_num_threads);
  stop_timer("1,4 interaction");

  DEBUG(6, "\tQMMM LJ exceptions");
  start_timer("LJ exceptions");
  m_outerloop.lj_exception_outerloop(topo, conf, sim, m_storage, m_rank, m_num_threads);
  stop_timer("LJ exceptions");
  
  // add long-range force

  return 0;
}

int interaction::QMMM_Nonbonded_Set::update_configuration(
                            const topology::Topology& topo,
                            configuration::Configuration& conf,
                            const simulation::Simulation& sim)
{
  const int ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & e = conf.current().energies;
  
  // use the IMPULSE method for multiple time stepping
  const unsigned int num_atoms = topo.num_atoms();
  math::VArray & force = conf.current().force;
  if (sim.param().multistep.steps > 1){
    int steps = sim.param().multistep.steps;
    if (sim.param().multistep.boost == 0)
      steps = 1;
    
    // only add when calculated
    if ((sim.steps() % steps) == 0){
      for(unsigned int i=0; i<num_atoms; ++i)
	      force(i) += steps * m_storage.force(i);
    }
    
  } else { // no multistep
    for(unsigned int i=0; i<num_atoms; ++i)
      force(i) += m_storage.force(i);
  }
  
  // (MULTISTEP: and keep energy constant)
  for (int i = 0; i < ljs; ++i) {
    for (int j = 0; j < ljs; ++j) {
      e.lj_energy[i][j] +=
              m_storage.energies.lj_energy[i][j];
      if (sim.param().force.force_groups) {
        for(unsigned int k = 0; k < num_atoms; ++k) {
          conf.special().force_groups[i][j][k] +=
                  m_storage.force_groups[i][j][k];
        }
      }
    }
  }
  if (sim.param().pcouple.virial) {
    DEBUG(7, "\tadd QMMM set virial");
  	conf.current().virial_tensor += m_storage.virial_tensor;
  }
  return 0;
}

int interaction::QMMM_Nonbonded_Set::init(const topology::Topology& topo
                                        , const configuration::Configuration& conf
                                        , const simulation::Simulation& sim
                                        , std::ostream& os
                                        , bool quiet)
{
  const unsigned num_atoms = topo.num_atoms();

  m_storage.force.resize(num_atoms);
  m_storage.energies.resize(
              unsigned(conf.current().energies.bond_energy.size())
            , unsigned(conf.current().energies.kinetic_energy.size()));
  
  if (sim.param().force.force_groups) {
    m_storage.force_groups.resize(unsigned(conf.current().energies.bond_energy.size()),
            std::vector<math::VArray>(unsigned(conf.current().energies.bond_energy.size()), 
            math::VArray(num_atoms, math::Vec(0.0, 0.0, 0.0))));
  }
  
  this->m_pairlist.resize(num_atoms);
  
  // guess the number of pairs
  const unsigned num_solute_atoms = topo.num_solute_atoms();
  unsigned num_qm_atoms = 0;
  for (unsigned i = 0; i < num_solute_atoms; ++i)
    num_qm_atoms += unsigned(topo.is_qm(i));


  const double vol = math::volume(conf.current().box, conf.boundary_type);
  unsigned pairs = 0;
  if (vol && sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    const double c3 = sim.param().qmmm.cutoff
                    * sim.param().qmmm.cutoff
                    * sim.param().qmmm.cutoff;
    
    pairs = num_qm_atoms * 
      unsigned(1.3 * num_atoms / vol * 4.0 / 3.0 * math::Pi * c3);
    
    if (pairs > num_atoms) pairs = num_atoms;

    if (!quiet)
      os << "\testimated pairlist size (QMMM) : "
	       << pairs << "\n";
  }
  else pairs = num_atoms;  
  this->m_pairlist.reserve(pairs);
  return 0;
}

int interaction::QMMM_Nonbonded_Set
::calculate_hessian(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    unsigned int atom_i, unsigned int atom_j,
		    math::Matrix & hessian){
  
  return m_outerloop.calculate_hessian(topo, conf, sim,
				       atom_i, atom_j, hessian,
				       m_pairlist);
}
