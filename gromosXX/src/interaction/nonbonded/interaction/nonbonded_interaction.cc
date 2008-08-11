/**
 * @file nonbonded_interaction.cc
 * template methods of Nonbonded_Interaction.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <simulation/parameter.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>

#include <interaction/nonbonded/interaction/storage.h>

#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>
#include <interaction/nonbonded/interaction/nonbonded_set_interface.h>
#include <interaction/nonbonded/interaction/nonbonded_set.h>
#include <interaction/nonbonded/interaction/vgrid_nonbonded_set.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_term.h>
#include <interaction/nonbonded/interaction/eds_nonbonded_term.h>

#include <interaction/nonbonded/interaction/perturbed_nonbonded_pair.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h>
#include <interaction/nonbonded/interaction/eds_nonbonded_outerloop.h>

#include <interaction/nonbonded/interaction/perturbed_nonbonded_set.h>
#include <interaction/nonbonded/interaction/eds_nonbonded_set.h>

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
    m_set_size(1)
{
  m_pairlist_algorithm->timer_pointer(&m_timer);
}

/**
 * Destructor.
 * @bug change destruction of nonbonded set to be standard - conform!
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

  m_timer.start();

  // check if we want to calculate nonbonded
  // might not be necessary if multiple time-stepping is enabled

  int steps = sim.param().multistep.steps;
  if (steps == 0) steps = 1;
  
  // std::cerr << "Nonbonded: steps = " << steps << std::endl;
  configuration::Configuration *exp_conf = NULL;
  if ((sim.steps() % steps) == 0){

    // std::cerr << "\tMULTISTEP: full non-bonded calculation" << std::endl;
    
    ////////////////////////////////////////////////////
    // multiple unit cell
    ////////////////////////////////////////////////////
    
    if (sim.param().multicell.multicell){
      
      DEBUG(6, "nonbonded: MULTICELL");
      exp_conf = new configuration::Configuration();
      expand_configuration(topo, conf, sim, *exp_conf);
      DEBUG(7, "\tmulticell conf: pos.size()=" << exp_conf->current().pos.size());

      // shared memory do this only once
      if(m_pairlist_algorithm->prepare(topo.multicell_topo(), *exp_conf, sim))
	return 1;
      
      // have to do all from here (probably it's only one,
      // but then maybe it's clearer like it is...)
      for(int i=0; i < m_set_size; ++i){
	m_nonbonded_set[i]->calculate_interactions(topo.multicell_topo(),
						   *exp_conf, 
						   sim);
      }
    }
    ////////////////////////////////////////////////////
    // end of MULTICELL
    ////////////////////////////////////////////////////
    else{
      
      // shared memory do this only once
      if (m_pairlist_algorithm->prepare(topo, conf, sim))
	return 1;
      
      // have to do all from here (probably it's only one,
      // but then maybe it's clearer like it is...)
      for(int i=0; i < m_set_size; ++i){
	m_nonbonded_set[i]->calculate_interactions(topo, conf, sim);
      }
    }

    ////////////////////////////////////////////////////
    // end of multiple time stepping: calculate
    ////////////////////////////////////////////////////
  }
  else{
    // std::cerr << "\tMULTISTEP: no non-bonded calculation" << std::endl;
  }
  
  DEBUG(6, "sets are done, adding things up...");
  store_set_data(topo, conf, sim);
  
  if (sim.param().multicell.multicell)  
    reduce_configuration(topo, conf, sim, *exp_conf);

  ////////////////////////////////////////////////////
  // printing pairlist
  ////////////////////////////////////////////////////
  if (sim.param().pairlist.print &&
      (!(sim.steps() % sim.param().pairlist.skip_step))){
    DEBUG(7, "print pairlist...");
    std::cerr << "printing pairlist!" << std::endl;
    if (sim.param().pairlist.grid != 2)
      print_pairlist(topo, conf, sim);
  }
  
  DEBUG(6, "Nonbonded_Interaction::calculate_interactions done");
  m_timer.stop();
  
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
int interaction::Nonbonded_Interaction::calculate_hessian
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 unsigned int atom_i, unsigned int atom_j,
 math::Matrix & hessian
 )
{
  std::vector<Nonbonded_Set_Interface *>::iterator
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
					     std::ostream & os,
					     bool quiet)
{ 
  if (!quiet)
    os << "Nonbonded interaction\n";

  // initialise the pairlist...
  m_pairlist_algorithm->init(topo, conf, sim, os, quiet);
  
  if (sim.param().nonbonded.method != simulation::el_reaction_field) {
    conf.lattice_sum().init(topo, sim);
  }

  DEBUG(15, "nonbonded_interaction::initialize");
  m_nonbonded_set.clear();
  
  if (sim.param().perturbation.perturbation){
    
    // std::cerr << "creating " << m_set_size << " Perturbed_Nonbonded_Sets" << std::endl;

    for(int i=0; i<m_set_size; ++i)
      m_nonbonded_set.push_back(new Perturbed_Nonbonded_Set(*m_pairlist_algorithm,
							    m_parameter, i, m_set_size));
  }
  else if (sim.param().eds.eds){
    DEBUG(16, "creating EDS nonbonded set");
    for(int i=0; i<m_set_size; ++i){
      m_nonbonded_set.push_back(new Eds_Nonbonded_Set(*m_pairlist_algorithm,
                                                      m_parameter, i, m_set_size));
      DEBUG(16, "pushed back EDS nonbonded set");
    }
      
  }
  else{
    for(int i=0; i<m_set_size; ++i){
      if (sim.param().pairlist.grid == 2){
	// std::cerr << "creating VGrid_Nonbonded_Set" << std::endl;
	m_nonbonded_set.push_back(new VGrid_Nonbonded_Set(*m_pairlist_algorithm, 
							  m_parameter, i, m_set_size));
      }
      else{
	// std::cerr << "creating Nonbonded_Set" << std::endl;
	m_nonbonded_set.push_back(new Nonbonded_Set(*m_pairlist_algorithm, 
						    m_parameter, i, m_set_size));
      }
    }
  }

  std::vector<Nonbonded_Set_Interface *>::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();

  bool q = quiet;
  for( ; it != to; ++it){
    if (sim.param().multicell.multicell)
      (*it)->init(topo.multicell_topo(), conf, sim, os, quiet);
    else
      (*it)->init(topo, conf, sim, os, q);
    // only print first time...
    q = true;
  }

  if (!quiet)
    os << "\n";
  
  if (check_spc_loop(topo, conf, sim, os, quiet) != 0) {
    io::messages.add("SPC loop was disabled", "Nonbonded_Interaction",
            io::message::warning);
  }
  DEBUG(9, "nonbonded init done");
  return 0;
}


//***************************************************************************
// helper functions 
//***************************************************************************

int interaction::Nonbonded_Interaction::check_spc_loop
(
 topology::Topology const & topo,
 configuration::Configuration const & conf,
 simulation::Simulation & sim,
 std::ostream & os,
 bool quiet)
{
  DEBUG(7, "checking for spc interaction loops");
  DEBUG(10, " param: special_loop = " << sim.param().force.special_loop);
  
  if (sim.param().force.special_loop ==  simulation::special_loop_off){
    DEBUG(8, "standard loops, user request");
    // sim.param().force.special_loop = 0;
    if (!quiet)
      os << "\tusing standard solvent loops (user request)\n";
    return 0;
  }
  
  if (sim.param().pairlist.atomic_cutoff){
    sim.param().force.special_loop =  simulation::special_loop_off;
    if (!quiet)
      os << "\tusing standard solvent loops (atomic cutoff)\n";
    return 1;
  }
  
  DEBUG(10, "num_solvents = " << topo.num_solvents());
  DEBUG(10, "molecules = " << topo.num_solvent_molecules(0));
  DEBUG(10, "atoms = " << topo.num_solvent_atoms(0));
  
  if (topo.num_solvents() != 1 ||
      topo.num_solvent_molecules(0) < 1 ||
      topo.num_solvent_atoms(0) / topo.num_solvent_molecules(0) != 3
      ){
  
    DEBUG(10, "standard loops...");
    
    sim.param().force.special_loop =  simulation::special_loop_off;
    if (!quiet)
      if (topo.num_solvents() > 0){
	if (topo.num_solvent_molecules(0) == 0){
	  os << "\tusing standard solvent loops (no solvent present!)\n\n";
	}
	else{
	  os << "\tusing standard solvent loops (num solvents doesn't match)\n"
	     << "\t\tnum solvents: " 
	     << topo.num_solvents() << "\n"
	     << "\t\tsolvent atoms: " 
	     << topo.num_solvent_atoms(0) / topo.num_solvent_molecules(0) << "\n"
	     << "\t\tmolecules: " << topo.num_solvent_molecules(0) << "\n\n";
	}
      }
      else{
	os << "\tusing standard solvent loops (no solvent in topology!)\n\n";
      }
    return 1;
  }
  
  DEBUG(10, "checking charges...");

  // check charges
  if (topo.charge()(topo.num_solute_atoms()) != -0.82 ||
      topo.charge()(topo.num_solute_atoms()+1) != 0.41 ||
      topo.charge()(topo.num_solute_atoms()+2) != 0.41){
    
    DEBUG(10, "charges don't match, standard loops");
    sim.param().force.special_loop =  simulation::special_loop_off;
    if (!quiet)
	os << "\tusing standard solvent loops (charges don't match)\n"
		  << "\t\tO  : " << topo.charge()(topo.num_solute_atoms()) << "\n"
		  << "\t\tH1 : " << topo.charge()(topo.num_solute_atoms()+1) << "\n"
		  << "\t\tH2 : " << topo.charge()(topo.num_solute_atoms()+2) << "\n\n";
	  
    return 1;
  }
  
  // check lj parameters
  DEBUG(10, "checking LJ parameter...");
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
    
    DEBUG(10, "don't match, force standard loop");
    sim.param().force.special_loop = simulation::special_loop_off;
    if (!quiet)
      os << "\tusing standard solvent loops (van der Waals parameter don't match)\n";
    return 1;
  }
  
  // check four_pi_eps_i
  DEBUG(10, "checking (4 pi eps0)^-1 ...");
  if(math::four_pi_eps_i!=138.9354){
    DEBUG(10, " does not match, force standard loop");
    sim.param().force.special_loop = simulation::special_loop_off;
    if(!quiet)
      os << "\tusing standard solvent loops ((4 pi eps0)^-1 does not match)\n";
    return 1;
  }
   
  // check energy groups
  DEBUG(10, "checking energy groups ...");

  if(topo.atom_energy_group(topo.num_solute_atoms())!=
          topo.atom_energy_group(topo.num_atoms()-1)){
    DEBUG(10, "- incompatible.");
    sim.param().force.special_loop = simulation::special_loop_off;
    if(!quiet)
      os << "\tusing standard solvent loops (energy group partitioning incompatible).\n"
              << "\tAll solvent atoms must be in one single energy group. ";
    return 1;
  }
    
  DEBUG(10, "happy to force spc loops");
  sim.param().force.special_loop = simulation::special_loop_off;
  if (!quiet)
    os << "\tusing spc solvent loops\n";
  
  return 0;
    
}

/**
 * store data from sets into the configuration
 */
void interaction::Nonbonded_Interaction::store_set_data
(
 topology::Topology const & topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim
 )
{
  std::vector<Nonbonded_Set_Interface *>::iterator
    it = m_nonbonded_set.begin(),
    to = m_nonbonded_set.end();
  
  // add the forces, energies, virial...
  for( ; it != to; ++it){
    DEBUG(7, "adding forces from set " << it - m_nonbonded_set.begin());
    (*it)->update_configuration(topo, conf, sim);
  }
}

/**
 * expand a configuration for
 * multiple unit cell
 * simulations
 */
void interaction::Nonbonded_Interaction::expand_configuration
(
 topology::Topology const & topo,
 configuration::Configuration const & conf,
 simulation::Simulation & sim,
 configuration::Configuration & exp_conf
 )
{

  DEBUG(7, "expanding configuration");
  
  const int mul = sim.param().multicell.x * sim.param().multicell.y * sim.param().multicell.z;
  DEBUG(8, "\tmul=" << mul);
  
  assert(topo.multicell_topo().num_atoms() == topo.num_atoms() * mul);
  // resize the configuration
  // exp_conf.resize(topo.multicell_topo().num_atoms());

  exp_conf.boundary_type = conf.boundary_type;
  exp_conf.current().box(0) = sim.param().multicell.x * conf.current().box(0);
  exp_conf.current().box(1) = sim.param().multicell.y * conf.current().box(1);
  exp_conf.current().box(2) = sim.param().multicell.z * conf.current().box(2);
  
  exp_conf.init(topo.multicell_topo(), sim.param(), false);
  DEBUG(10, "\texp_conf initialised");
  
  math::Vec shift(0.0);
  int exp_i=0;

  // SOLUTE
  DEBUG(10, "\tsolute");
  for(int z=0; z<sim.param().multicell.z; ++z){
    for(int y=0; y<sim.param().multicell.y; ++y){
      for(int x=0; x<sim.param().multicell.x; ++x){

	shift = x * conf.current().box(0) + 
	  y * conf.current().box(1) + 
	  z * conf.current().box(2);

	// this should be the NORMAL topo!
	for(unsigned int i=0; i<topo.num_solute_atoms(); ++i, ++exp_i){
	  
	  assert(exp_conf.current().pos.size() > unsigned(exp_i));
	  assert(conf.current().pos.size() > i);
	  
	  exp_conf.current().pos(exp_i) = conf.current().pos(i) + shift;
          exp_conf.current().posV(exp_i) = conf.current().posV(i);
	  // exp_conf.old().pos(exp_i) = conf.old().pos(i) + shift;
	}
      }
    }
  }

  // SOLVENT
  for(int z=0; z<sim.param().multicell.z; ++z){
    for(int y=0; y<sim.param().multicell.y; ++y){
      for(int x=0; x<sim.param().multicell.x; ++x){

	shift = x * conf.current().box(0) + 
	  y * conf.current().box(1) + 
	  z * conf.current().box(2);

	for(unsigned int i=topo.num_solute_atoms(); i<topo.num_atoms(); ++i, ++exp_i){
	  
	  assert(exp_conf.current().pos.size() > unsigned(exp_i));
	  assert(conf.current().pos.size() > i);
	  
	  exp_conf.current().pos(exp_i) = conf.current().pos(i) + shift;
          exp_conf.current().posV(exp_i) = conf.current().posV(i);
	  // exp_conf.old().pos(exp_i) = conf.old().pos(i) + shift;
	}
      }
    }
  }
  
  exp_conf.current().energies.zero();
  exp_conf.current().perturbed_energy_derivatives.zero();
  exp_conf.current().virial_tensor = 0.0;

}

/**
 * reduce a configuration for
 * multiple unit cell
 * simulations
 */
void interaction::Nonbonded_Interaction::reduce_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 configuration::Configuration & exp_conf
 )
 {
  // add one-four, rf excluded etc... all those things that go directly into 
  // the configuration and not into storages of the sets
  
  // reduce the forces
  const unsigned int cells = (sim.param().multicell.x * sim.param().multicell.y * sim.param().multicell.z);
  unsigned int i = 0;
  for(; i < topo.num_solute_atoms(); ++i) {
    conf.current().force(i) += exp_conf.current().force(i);
    conf.current().posV(i) += exp_conf.current().posV(i);
  }
  
  // one cell is is already contained in i!! -> cells - 1
  const unsigned int offset = topo.num_solute_atoms() * (cells - 1);
  for(; i < topo.num_atoms(); ++i) {
    conf.current().force(i) += exp_conf.current().force(offset + i);
    conf.current().posV(i) += exp_conf.current().posV(offset + i);
  }
  
  const int ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & e = conf.current().energies;
  configuration::Energy & exp_e = exp_conf.current().energies;
  configuration::Energy & pe = conf.current().perturbed_energy_derivatives;
  configuration::Energy & exp_pe = exp_conf.current().perturbed_energy_derivatives;

  const double cells_i = 1.0 / cells;
  
  // reduce the energies
  for(int i = 0; i < ljs; ++i){
    for(int j = 0; j < ljs; ++j){
      e.lj_energy[i][j] += exp_e.lj_energy[i][j] * cells_i;
      e.crf_energy[i][j] += exp_e.crf_energy[i][j] * cells_i;
      pe.lj_energy[i][j] += exp_pe.lj_energy[i][j] * cells_i;
      pe.crf_energy[i][j] += exp_pe.crf_energy[i][j] * cells_i;
    }
  }
  // reduce the virial
  if (sim.param().pcouple.virial){
    DEBUG(7, "\tadd set virial");
    
    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){
        
        conf.current().virial_tensor(i, j) +=
        exp_conf.current().virial_tensor(i, j) * cells_i;
      }
    }
  }
}


/**
 * print the pairlist
 */
int interaction::Nonbonded_Interaction::print_pairlist
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::ostream & os
 )
{
  DEBUG(4, "Nonbonded_Interaction::print_pairlist");
  
  Pairlist temp_solute, temp_solvent;
  temp_solute.resize(topo.num_atoms());
  temp_solvent.resize(topo.num_atoms());
  
  for(unsigned int atom_i = 0; atom_i < topo.num_atoms(); ++atom_i){
    
    for(int i=0; i < m_set_size; ++i){

      assert (m_nonbonded_set.size() > unsigned(i));
      assert (m_nonbonded_set[i]->pairlist().solute_short.size() > atom_i);
      assert (m_nonbonded_set[i]->pairlist().solvent_short.size() > atom_i);
      
      for(unsigned int atom_j = 0;
	  atom_j < m_nonbonded_set[i]->pairlist().solute_short[atom_i].size();
	  ++atom_j){
	
	assert(temp_solute.size() > atom_i);
	assert(temp_solute.size() > m_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]);

	if (m_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j] < atom_i)
	  temp_solute[m_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]].push_back(atom_i);
	else
	  temp_solute[atom_i].push_back(m_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]);
      }
      for(unsigned int atom_j = 0;
	  atom_j < m_nonbonded_set[i]->pairlist().solvent_short[atom_i].size();
	  ++atom_j){
	
	assert(temp_solvent.size() > atom_i);
	assert(temp_solvent.size() > m_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]);

	if (m_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j] < atom_i)
	  temp_solvent[m_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]].push_back(atom_i);
	else
	  temp_solvent[atom_i].push_back(m_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]);
      }
    }
  }
  
  os << temp_solute << std::endl
     << temp_solvent << std::endl;
  return 0;
}
