/**
 * @file nonbonded_set.cc
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

#include "../../../interaction/nonbonded/interaction/latticesum.h"

#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"

#include "../../../math/periodicity.h"
#include "../../../math/volume.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_set.h"

#include "../../../util/debug.h"

#include "../../../configuration/energy.h"

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
  const bool pairlist_update = !(sim.steps() % sim.param().pairlist.skip_step);
  if(pairlist_update){
    DEBUG(6, "\tdoing longrange...");
    
    //====================
    // create a pairlist
    //====================
    
    // zero the longrange forces, energies, virial
    m_longrange_storage.zero();
    m_pairlist_alg.update(topo, conf, sim, 
			  pairlist(),
			  m_rank, topo.num_atoms(), m_num_threads);
  }

  // Print pairlist
  //if (sim.steps() ==5){
  //  std::cout << pairlist().solute_long << std::endl;
  //  std::cout << pairlist().solute_short << std::endl;
  //  std::cout << pairlist().solvent_long << std::endl;
  //  std::cout << pairlist().solvent_short << std::endl;
  //}

  if (sim.param().polarise.cos) {
    //===============
    // polarisation
    //===============

    // calculate explicit polarisation of the molecules
    DEBUG(6, "\texplicit polarisation");
    start_subtimer("explicit polarisation");
    
    m_outerloop.electric_field_outerloop(topo, conf, sim, m_pairlist,
				       m_storage, m_longrange_storage, m_rank);
    stop_subtimer("explicit polarisation");
  }
  
  //DEBUG(1, "field 1 [" << m_storage.electric_field(0)(0) << " " << m_storage.electric_field(0)(1) << " "  << m_storage.electric_field(0)(2) << "]");
  //DEBUG(1, "lr field 1 [" << m_longrange_storage.electric_field(0)(0) << " " << m_longrange_storage.electric_field(0)(1) << " "  << m_longrange_storage.electric_field(0)(2) << "]");
  
  if (pairlist_update) {
    DEBUG(8, "doing longrange calculation");
    //start_subtimer("longrange interactions (total)");    //This timer counts all nonbonded long range interactions
    // in case of LS only LJ are calculated

    m_outerloop.lj_crf_outerloop(topo, conf, sim,
            m_pairlist.solute_long, m_pairlist.solvent_long,
            m_longrange_storage, true /*longrange!*/, m_pairlist_alg.timer(),
            m_rank == 0);
    //stop_subtimer("longrange interactions (total)");
  }


  // calculate forces / energies
  DEBUG(6, "\tshort range interactions");
  //start_subtimer("shortrange interactions (total)");

  if (sim.param().force.interaction_function == simulation::lj_ls_func) {
    start_subtimer("shortrange ls_real_outerloop");
    m_outerloop.ls_real_outerloop(topo, conf, sim,
            m_pairlist.solute_short, m_pairlist.solvent_short,
            m_storage, m_rank, m_num_threads);
    stop_subtimer("shortrange ls_real_outerloop");
  } else {
    m_outerloop.lj_crf_outerloop(topo, conf, sim,
            m_pairlist.solute_short, m_pairlist.solvent_short,
            m_storage, false, m_pairlist_alg.timer(),
            m_rank == 0);
  }

  //stop_subtimer("shortrange interactions (total)");

  // calculate sasa interaction
  if (sim.param().sasa.switch_sasa) {
    DEBUG(6, "\tsasa energy");
    start_subtimer("SASA energy");
    m_outerloop.sasa_outerloop(topo, conf, sim, m_storage);
    stop_subtimer("SASA energy");
  }

  // calculate k-space energy and forces
  switch (sim.param().nonbonded.method) {
    case simulation::el_ewald :
    {
      start_subtimer("k-space Ewald");
      DEBUG(6, "\tlong range electrostatics: Ewald");
      if (sim.param().pcouple.scale != math::pcouple_off || sim.steps() == 0) {
        // the box may have changed. Recalculate the k space
        configuration::calculate_k_space(topo, conf, sim, m_rank, m_num_threads);
      }
      // do the longrange calculation in k space
      m_outerloop.ls_ewald_kspace_outerloop(topo, conf, sim, m_storage,
              m_rank, m_num_threads);
      stop_subtimer("k-space Ewald");
      break;
    }
    case simulation::el_p3m :
    {
      start_subtimer("k-space P3M");
      DEBUG(6, "\tlong range electrostatics: P3M");
      if (sim.param().pcouple.scale != math::pcouple_off || sim.steps() == 0) {
        // check whether we have to recalculate the influence function
      }
      // do the longrange calculation in k space
      bool is_ok = true;
      m_outerloop.ls_p3m_kspace_outerloop(topo, conf, sim, m_storage,
              m_rank, m_num_threads, m_pairlist_alg.timer(), is_ok);
      stop_subtimer("k-space P3M");
  // !!! crash if not ok
        if ( is_ok == false ) {
          return 1;
        }
      break;
    }
    default:
    {
    } // doing reaction field
  }

  // calculate lattice sum self energy and A term
  // this has to be done after the k-space energy is calculated
  // as this calculation will also deliver a methodology
  // dependent A~_2 term if requested
  
  // check whether we have to calculate it at all in this step
  const bool calculate_lattice_sum_corrections =
          sim.param().pcouple.scale != math::pcouple_off || // NPT - every step
          !sim.steps() || // at the beginning of the simulation
          (sim.param().write.energy && sim.steps() % abs(sim.param().write.energy)) == 0; // energy output req.

  switch (sim.param().nonbonded.method) {
    case simulation::el_ewald :
    case simulation::el_p3m :
    {
      if (calculate_lattice_sum_corrections) {
        start_subtimer("LS self energy and A term");
        m_outerloop.ls_self_outerloop(topo, conf, sim, m_storage,
                m_rank, m_num_threads);
        stop_subtimer("LS self energy and A term");
      }
    }
    default:; // doing reaction field
  }

  DEBUG(6, "\t1,4 - interactions");
  start_subtimer("1,4 interaction");
  m_outerloop.one_four_outerloop(topo, conf, sim, m_storage,  m_rank, m_num_threads);
  stop_subtimer("1,4 interaction");

  start_subtimer("LJ exceptions");
  m_outerloop.lj_exception_outerloop(topo, conf, sim, m_storage,  m_rank, m_num_threads);
  stop_subtimer("LJ exceptions");

  // possibly do the RF contributions due to excluded atoms
  if (sim.param().nonbonded.rf_excluded) {
    DEBUG(7, "\tRF excluded interactions and self term");
    start_subtimer("RF excluded interaction");
    m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_storage, m_rank, m_num_threads);
    stop_subtimer("RF excluded interaction");
  }

  // single processor algorithms 
  if (m_rank == 0) {
    // calculate lattice sum surface energy and force
    
    switch (sim.param().nonbonded.method) {
      case simulation::el_ewald :
      case simulation::el_p3m :
      {
        start_subtimer("LS surface energy");
        m_outerloop.ls_surface_outerloop(topo, conf, sim, m_storage,
                m_rank, m_num_threads);
        stop_subtimer("LS surface energy");
      }
      default:
      {
      } // doing reaction field 
    }

    if (sim.param().polarise.cos) {
      DEBUG(6, "\tself-energy of polarisation");
      start_subtimer("polarisation self-energy");
      m_outerloop.self_energy_outerloop(topo, conf, sim, m_storage);
      stop_subtimer("polarisation self-energy");
    }
  }

  // add long-range force
  DEBUG(6, "\t(set) add long range forces");

  start_subtimer("longrange addition");
  m_storage.force += m_longrange_storage.force;
  
  // and long-range energies
  DEBUG(6, "\t(set) add long range energies");
  const unsigned int lj_e_size = unsigned(m_storage.energies.lj_energy.size());
  const unsigned int num_atoms = topo.num_atoms();

  for (unsigned int i = 0; i < lj_e_size; ++i) {
    for (unsigned int j = 0; j < lj_e_size; ++j) {
      m_storage.energies.lj_energy[i][j] +=
              m_longrange_storage.energies.lj_energy[i][j];
      m_storage.energies.crf_energy[i][j] +=
              m_longrange_storage.energies.crf_energy[i][j];
      m_storage.energies.shift_extra_orig[i][j] +=
              m_longrange_storage.energies.shift_extra_orig[i][j];
      m_storage.energies.shift_extra_phys[i][j] +=
              m_longrange_storage.energies.shift_extra_phys[i][j];
      if (sim.param().force.force_groups) {
        for (unsigned int k = 0; k < num_atoms; ++k) {
          m_storage.force_groups[i][j][k] +=
                  m_longrange_storage.force_groups[i][j][k];
        } // atoms
      } // if force groups
    }
  }
  //ORIOL_GAMD check
  if(sim.param().gamd.gamd){
    const unsigned int nr_agroups = unsigned(sim.param().gamd.igroups);
    for (unsigned int igroup = 0; igroup < nr_agroups; ++igroup) {
      m_storage.energies.gamd_potential_total[igroup] += m_longrange_storage.energies.gamd_potential_total[igroup];
      DEBUG(15, "Long range GAMD energys " << m_longrange_storage.energies.gamd_potential_total[igroup]);
      m_storage.force_gamd[igroup] += m_longrange_storage.force_gamd[igroup];
    }
  }

  // add longrange virial
  if (sim.param().pcouple.virial){
    DEBUG(6, "\t(set) add long range virial");
	m_storage.virial_tensor += m_longrange_storage.virial_tensor;
  //ORIOL_GAMD
  if(sim.param().gamd.gamd){
    const unsigned int nr_agroups = unsigned(sim.param().gamd.igroups);
    for (unsigned int igroup = 0; igroup < nr_agroups; ++igroup) {
      m_storage.virial_tensor_gamd[igroup] += m_longrange_storage.virial_tensor_gamd[igroup];
    }
  }
  }
  stop_subtimer("longrange addition");
  DEBUG(6, "\t(set) DONE");
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
  const unsigned int num_atoms = topo.num_atoms();
  math::VArray & f = conf.current().force;
  if (sim.param().multistep.steps > 1){
    int steps = sim.param().multistep.steps;
    if (sim.param().multistep.boost == 0)
      steps = 1;
    
    // only add when calculated
    if ((sim.steps() % steps) == 0){

      // std::cerr << "\tadding boosted (" << steps << ") non-bonded forces" << std::endl;

      for(unsigned int i=0; i<num_atoms; ++i)
	f(i) += steps * m_storage.force(i);
    }
    else{
      // std::cerr << "\tnot adding non-bonded forces" << std::endl;
    }
    
  } else{ // no multistep
    for(unsigned int i=0; i<num_atoms; ++i){
      f(i) += m_storage.force(i);
    }
      // ORIOL_GAMD
    if (sim.param().gamd.gamd){
      const unsigned int nr_agroups = unsigned(sim.param().gamd.igroups);
      for (unsigned int igroup=0; igroup < nr_agroups; ++igroup){
        DEBUG(5, "adding energy " << m_storage.energies.gamd_potential_total[igroup]);
        e.gamd_potential_total[igroup] += m_storage.energies.gamd_potential_total[igroup];
          for(unsigned int atomn=0; atomn<num_atoms; ++atomn){
            conf.special().gamd.total_force[igroup](atomn) +=  m_storage.force_gamd[igroup](atomn);
          } // loop over atoms
      } // loop over acceleration groups
    //energies
    } // end if
    DEBUG(5, "FINISHING UPATING FORCES");
  }

  
  // (MULTISTEP: and keep energy constant)
  for (int i = 0; i < ljs; ++i) {
    for (int j = 0; j < ljs; ++j) {
      e.lj_energy[i][j] +=
              m_storage.energies.lj_energy[i][j];
      e.crf_energy[i][j] +=
              m_storage.energies.crf_energy[i][j];
      e.shift_extra_orig[i][j] +=
              m_storage.energies.shift_extra_orig[i][j];
      e.shift_extra_phys[i][j] +=
              m_storage.energies.shift_extra_phys[i][j];
      e.ls_real_energy[i][j] +=
              m_storage.energies.ls_real_energy[i][j];
      e.ls_k_energy[i][j] +=
              m_storage.energies.ls_k_energy[i][j];
      if (sim.param().force.force_groups) {
        for(unsigned int k = 0; k < num_atoms; ++k) {
          conf.special().force_groups[i][j][k] +=
                  m_storage.force_groups[i][j][k];
        }
      }
    }
    e.self_energy[i] += m_storage.energies.self_energy[i];

  }


  // ANITA
  const unsigned int nr_lambdas = unsigned(m_storage.energies.A_lj_total.size());
  for(unsigned int i=0; i < nr_lambdas; ++i) {
    for(int j=0; j < ljs; ++j) {
      for(int k=0; k < ljs; ++k) {
        e.A_lj_energy[i][j][k] += 
            m_storage.energies.A_lj_energy[i][j][k];
        e.B_lj_energy[i][j][k] += 
            m_storage.energies.B_lj_energy[i][j][k];
        e.A_crf_energy[i][j][k] += 
            m_storage.energies.A_crf_energy[i][j][k];
        e.B_crf_energy[i][j][k] += 
            m_storage.energies.B_crf_energy[i][j][k];
      }
    }
  } //ANITA


  // no components in lattice sum methods!
  
  // lattice sum energies
  e.ls_kspace_total += m_storage.energies.ls_kspace_total;
  e.ls_self_total += m_storage.energies.ls_self_total;
  e.ls_a_term_total += m_storage.energies.ls_a_term_total;
  e.ls_surface_total += m_storage.energies.ls_surface_total;


  // (MULTISTEP: and the virial???)
  if (sim.param().pcouple.virial){
    DEBUG(7, "\tadd set virial");
  	conf.current().virial_tensor += m_storage.virial_tensor;
      // ORIOL_GAMD
      DEBUG(5, "UPDATE VIRIAL");
      if (sim.param().gamd.gamd){
        const unsigned int nr_agroups = unsigned(sim.param().gamd.igroups);
        for (unsigned int igroup=0; igroup < nr_agroups; ++igroup){
            for (int a = 0; a < 3; ++a) {
               for (int b = 0; b < 3; ++b) {
                 conf.special().gamd.virial_tensor[igroup](b, a) +=  m_storage.virial_tensor_gamd[igroup](b, a);
                 DEBUG(5, "GROUP is " << igroup);
                 DEBUG(5, conf.special().gamd.virial_tensor[igroup](b, a));
                 DEBUG(5, conf.special().gamd.virial_tensor_dihe[igroup](b, a));
                 DEBUG(5, conf.current().virial_tensor(b,a));
               }
            }
            DEBUG(5, "VIRIAL CICLED COMPLETED");
        } // loop over acceleration groups
      } // end if
  } // end virial if
  
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
	   unsigned(conf.current().energies.kinetic_energy.size()),
           unsigned(sim.param().precalclam.nr_lambdas)); //ANITA
  m_longrange_storage.energies.
    resize(unsigned(conf.current().energies.bond_energy.size()),
	   unsigned(conf.current().energies.kinetic_energy.size()),
           unsigned(sim.param().precalclam.nr_lambdas)); //ANITA

  // ORIOL_GAMD
  m_storage.force_gamd.resize(sim.param().gamd.igroups);
  m_storage.virial_tensor_gamd.resize(sim.param().gamd.igroups);
  m_storage.energies.gamd_potential_total.resize(sim.param().gamd.igroups);

  m_longrange_storage.force_gamd.resize(sim.param().gamd.igroups);
  m_longrange_storage.virial_tensor_gamd.resize(sim.param().gamd.igroups);
  m_longrange_storage.energies.gamd_potential_total.resize(sim.param().gamd.igroups);

  for(unsigned int i = 0; i < m_storage.force_gamd.size(); i++){
      m_storage.force_gamd[i].resize(topo.num_atoms());
      m_longrange_storage.force_gamd[i].resize(topo.num_atoms());
  }
  
  if (sim.param().force.force_groups) {
    m_storage.force_groups.resize(unsigned(conf.current().energies.bond_energy.size()),
            std::vector<math::VArray>(unsigned(conf.current().energies.bond_energy.size()), 
            math::VArray(num_atoms, math::Vec(0.0, 0.0, 0.0))));
    m_longrange_storage.force_groups.resize(unsigned(conf.current().energies.bond_energy.size()),
            std::vector<math::VArray>(unsigned(conf.current().energies.bond_energy.size()), 
            math::VArray(num_atoms, math::Vec(0.0, 0.0, 0.0))));
  }
  
  m_storage.electric_field.resize(num_atoms);
  m_longrange_storage.electric_field.resize(num_atoms);
  
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
      os << "\testimated pairlist size (per atom) : "
	 << pairs << "\n";
    
    pairlist().reserve(pairs);
  }
  
  return 0;
}

