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

#include <interaction/nonbonded/interaction/latticesum.h>

#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>

#include <math/periodicity.h>
#include <math/volume.h>

#include <interaction/nonbonded/interaction/nonbonded_set.h>

#include <util/debug.h>

#include <configuration/energy.h>

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

    // parallelisation using STRIDE:
    // chargegroup based pairlist can only use this one!!!!
    // TODO:
    // move decision to pairlist!!!
    m_pairlist_alg.update(topo, conf, sim, 
			  pairlist(),
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

  if (sim.param().polarize.cos) {
    //===============
    // polarization
    //===============

    // calculate explicit polarization of the molecules
    DEBUG(6, "\texplicit polarization");
    if (m_rank == 0)
      m_pairlist_alg.timer().start("explicit polarization");
    
    m_outerloop.electric_field_outerloop(topo, conf, sim, m_pairlist,
				       m_storage, m_longrange_storage, m_rank);
    if (m_rank == 0)
      m_pairlist_alg.timer().stop("explicit polarization");
  }
  
  //DEBUG(1, "field 1 [" << m_storage.electric_field(0)(0) << " " << m_storage.electric_field(0)(1) << " "  << m_storage.electric_field(0)(2) << "]");
  //DEBUG(1, "lr field 1 [" << m_longrange_storage.electric_field(0)(0) << " " << m_longrange_storage.electric_field(0)(1) << " "  << m_longrange_storage.electric_field(0)(2) << "]");
  
  if (pairlist_update) {
    DEBUG(8, "doing longrange calculation");
    if (m_rank == 0)
      m_pairlist_alg.timer().start("longrange");
    // in case of LS only LJ are calculated

    m_outerloop.lj_crf_outerloop(topo, conf, sim,
            m_pairlist.solute_long, m_pairlist.solvent_long,
            m_longrange_storage);
    if (m_rank == 0)
      m_pairlist_alg.timer().stop("longrange");
  }


  // calculate forces / energies
  DEBUG(6, "\tshort range interactions");
  if (m_rank == 0)
    m_pairlist_alg.timer().start("shortrange");

  if (sim.param().force.interaction_function == simulation::lj_ls_func) {
    m_outerloop.ls_real_outerloop(topo, conf, sim,
            m_pairlist.solute_short, m_pairlist.solvent_short,
            m_storage, m_rank, m_num_threads);
  } else {
    m_outerloop.lj_crf_outerloop(topo, conf, sim,
            m_pairlist.solute_short, m_pairlist.solvent_short,
            m_storage);
  }

  if (m_rank == 0)
    m_pairlist_alg.timer().stop("shortrange");

  // calculate k-space energy and forces
  if (m_rank == 0)
    m_pairlist_alg.timer().start("k-space");
  switch (sim.param().nonbonded.method) {
    case simulation::el_ewald :
    {
      DEBUG(6, "\tlong range electrostatics: Ewald");
      if (sim.param().pcouple.scale != math::pcouple_off || sim.steps() == 0) {
        // the box may have changed. Recalculate the k space
        configuration::calculate_k_space(topo, conf, sim);
      }
      // do the longrange calculation in k space
      m_outerloop.ls_ewald_kspace_outerloop(topo, conf, sim, m_storage,
              m_rank, m_num_threads);
      break;
    }
    case simulation::el_p3m :
    {
      DEBUG(6, "\tlong range electrostatics: P3M");
      if (sim.param().pcouple.scale != math::pcouple_off || sim.steps() == 0) {
        // check whether we have to recalculate the influence function
      }
      // do the longrange calculation in k space
      m_outerloop.ls_p3m_kspace_outerloop(topo, conf, sim, m_storage,
              m_rank, m_num_threads, m_pairlist_alg.timer());
      break;
    }
    default:
    {
    } // doing reaction field
  }
  if (m_rank == 0)
    m_pairlist_alg.timer().stop("k-space");
  
  // single processor algorithms (1-4,...)
  if (m_rank == 0) {
    // calculate lattice sum self energy and A term
    // this has to be done after the k-space energy is calculated
    // as this calculation will also deliver a methodology
    // dependent A~_2 term if requested

    m_pairlist_alg.timer().start("ls self energy and A term");
    switch (sim.param().nonbonded.method) {
      case simulation::el_ewald :
      case simulation::el_p3m :
      {
        //  NPT                                            ||  NVT
        if (sim.param().pcouple.scale != math::pcouple_off || sim.steps() == 0) {
          m_outerloop.ls_self_outerloop(topo, conf, sim, m_storage,
                  m_rank, m_num_threads);
        } else {
          // copy from previous step
          conf.current().energies.ls_self_total = conf.old().energies.ls_self_total;
          conf.current().energies.ls_a_term_total = conf.old().energies.ls_a_term_total;
        }
      }
      default:; // doing reaction field
    }

    m_pairlist_alg.timer().stop("ls self energy and A term");

    // calculate lattice sum surface energy and force

    m_pairlist_alg.timer().start("ls surface energy");
    switch (sim.param().nonbonded.method) {
      case simulation::el_ewald :
      case simulation::el_p3m :
      {
        m_outerloop.ls_surface_outerloop(topo, conf, sim, m_storage,
                m_rank, m_num_threads);
      }
      default:
      {
      } // doing reaction field 
    }

    m_pairlist_alg.timer().stop("ls surface energy");



    DEBUG(6, "\t1,4 - interactions");
    m_pairlist_alg.timer().start("1,4 interaction");
    m_outerloop.one_four_outerloop(topo, conf, sim, m_storage);
    m_pairlist_alg.timer().stop("1,4 interaction");

    // possibly do the RF contributions due to excluded atoms
    if (sim.param().nonbonded.rf_excluded) {
      DEBUG(7, "\tRF excluded interactions and self term");
      m_pairlist_alg.timer().start("RF excluded interaction");
      m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_storage);
      m_pairlist_alg.timer().stop("RF excluded interaction");
    }

    if (sim.param().polarize.cos) {
      DEBUG(6, "\tself-energy of polarization");
      m_pairlist_alg.timer().start("polarization self-energy");
      m_outerloop.self_energy_outerloop(topo, conf, sim, m_storage);
      m_pairlist_alg.timer().stop("polarization self-energy");
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

	m_storage.virial_tensor += m_longrange_storage.virial_tensor;
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

  const unsigned int cells = (sim.param().multicell.x * sim.param().multicell.y * sim.param().multicell.z);
  
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
    if (sim.param().multicell.multicell) {
      unsigned int i = 0;
      for(; i < topo.num_solute_atoms(); ++i)
        conf.current().force(i) += m_storage.force(i);
      
      // one cell is is already contained in i!! -> cells - 1
      const unsigned int offset = topo.num_solute_atoms() * (cells - 1); 
      for(; i < topo.num_atoms(); ++i) {
        conf.current().force(i) += m_storage.force(offset + i);
      }
      
    } else {
      for(unsigned int i=0; i<topo.num_atoms(); ++i)
        conf.current().force(i) += m_storage.force(i);
    }
  }
  
  // (MULTISTEP: and keep energy constant)
  const double cells_i = 1.0 / cells;
  for(int i = 0; i < ljs; ++i){
    for(int j = 0; j < ljs; ++j){
      
      e.lj_energy[i][j] += 
	m_storage.energies.lj_energy[i][j] * cells_i;
      e.crf_energy[i][j] += 
	m_storage.energies.crf_energy[i][j] * cells_i;
      e.ls_real_energy[i][j] += 
	m_storage.energies.ls_real_energy[i][j] * cells_i;
      e.ls_k_energy[i][j] += 
	m_storage.energies.ls_k_energy[i][j] * cells_i;
    }
    e.self_energy[i] +=  m_storage.energies.self_energy[i] * cells_i;
  }
  // no components in lattice sum methods!
  
  // lattice sum energies
  e.ls_kspace_total += m_storage.energies.ls_kspace_total;
  e.ls_self_total += m_storage.energies.ls_self_total;
  e.ls_a_term_total += m_storage.energies.ls_a_term_total;
  e.ls_surface_total += m_storage.energies.ls_surface_total;


  // (MULTISTEP: and the virial???)
  if (sim.param().pcouple.virial){
    DEBUG(7, "\tadd set virial");
  	conf.current().virial_tensor += m_storage.virial_tensor * cells_i;
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
      os << "\n\testimated pairlist size (per atom) : "
	 << pairs << "\n";
    
    pairlist().reserve(pairs);
  }
  
  return 0;
}

