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
#include <interaction/nonbonded/interaction/nonbonded_set_interface.h>
#include <interaction/nonbonded/interaction/nonbonded_set.h>
#include <interaction/nonbonded/interaction/vgrid_nonbonded_set.h>

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
    m_longrange_timing(0.0),
    m_set_size(1)
{
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

  assert((sim.param().force.spc_loop <= 0) || 
	 (!sim.param().pairlist.grid && !sim.param().pairlist.atomic_cutoff));

  const double nonbonded_start = util::now();

  // multiple unit cell
  if (sim.param().multicell.multicell){
    
    DEBUG(6, "nonbonded: MULTICELL");
    configuration::Configuration exp_conf;
    expand_configuration(topo, conf, sim, exp_conf);
    DEBUG(7, "\tmulticell conf: pos.size()=" << exp_conf.current().pos.size());
    
    // shared memory do this only once
    if(m_pairlist_algorithm->prepare(topo.multicell_topo(), exp_conf, sim))
      return 1;

    // have to do all from here (probably it's only one,
    // but then maybe it's clearer like it is...)
    for(int i=0; i < m_set_size; ++i){
      m_nonbonded_set[i]->calculate_interactions(topo.multicell_topo(),
						 exp_conf, 
						 sim);
    }
  }
  else{ // no MULTICELL
    
    // shared memory do this only once
    if (m_pairlist_algorithm->prepare(topo, conf, sim))
      return 1;
    
    // have to do all from here (probably it's only one,
    // but then maybe it's clearer like it is...)
    for(int i=0; i < m_set_size; ++i){
      m_nonbonded_set[i]->calculate_interactions(topo, conf, sim);
    }

    DEBUG(6, "sets are done, adding things up...");
    store_set_data(topo, conf, sim);
    
    DEBUG(7, "print pairlist...");
    if (sim.param().pairlist.print &&
	(!(sim.steps() % sim.param().pairlist.skip_step))){
      
      std::cerr << "printing pairlist!" << std::endl;
      if (sim.param().pairlist.grid != 2)
	print_pairlist(topo, conf, sim);
    }
  }
  
  DEBUG(6, "Nonbonded_Interaction::calculate_interactions done");
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

  DEBUG(15, "nonbonded_interaction::initialize");
  m_nonbonded_set.clear();
  
  if (sim.param().perturbation.perturbation){
    
    // std::cerr << "creating " << m_set_size << " Perturbed_Nonbonded_Sets" << std::endl;

    for(int i=0; i<m_set_size; ++i)
      m_nonbonded_set.push_back(new Perturbed_Nonbonded_Set(*m_pairlist_algorithm,
							    m_parameter, i, m_set_size));
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
  
  check_spc_loop(topo, conf, sim, os, quiet);
  DEBUG(9, "nonbonded init done");
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
 std::ostream & os,
 bool quiet)
{
  DEBUG(7, "checking for spc interaction loops");
  DEBUG(10, " param: spc_loop = " << sim.param().force.spc_loop);
  

  if (sim.param().force.spc_loop == -1){
    DEBUG(8, "standard loops, user request");
    // sim.param().force.spc_loop = 0;
    if (!quiet)
      os << "\tusing standard solvent loops (user request)\n";
    return;
  }

  if (sim.param().pairlist.grid){
    sim.param().force.spc_loop = 0;
    if (!quiet)
      os << "\tusing standard solvent loops (grid based pairlist)\n";
    return;
  }
  
  DEBUG(10, "num_solvents = " << topo.num_solvents());
  DEBUG(10, "molecules = " << topo.num_solvent_molecules(0));
  DEBUG(10, "atoms = " << topo.num_solvent_atoms(0));
  
  if (topo.num_solvents() != 1 ||
      topo.num_solvent_molecules(0) < 1 ||
      topo.num_solvent_atoms(0) / topo.num_solvent_molecules(0) != 3
      ){
  
    DEBUG(10, "standard loops...");
    
    sim.param().force.spc_loop = 0;
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
    return;
  }
  
  DEBUG(10, "checking charges...");

  // check charges
  if (topo.charge()(topo.num_solute_atoms()) != -0.82 ||
      topo.charge()(topo.num_solute_atoms()+1) != 0.41 ||
      topo.charge()(topo.num_solute_atoms()+2) != 0.41){
    
    DEBUG(10, "charges don't match, standard loops");
    sim.param().force.spc_loop = 0;
    if (!quiet)
	os << "\tusing standard solvent loops (charges don't match)\n"
		  << "\t\tO  : " << topo.charge()(topo.num_solute_atoms()) << "\n"
		  << "\t\tH1 : " << topo.charge()(topo.num_solute_atoms()+1) << "\n"
		  << "\t\tH2 : " << topo.charge()(topo.num_solute_atoms()+2) << "\n\n";
	  
    return;
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
    sim.param().force.spc_loop = 0;
    if (!quiet)
      os << "\tusing standard solvent loops (van der Waals parameter don't match)\n";
    return;
  }
  
  DEBUG(10, "happy to force spc loops");
  sim.param().force.spc_loop = 1;
  if (!quiet)
    os << "\tusing spc solvent loops\n";
    return;
    
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
	  // exp_conf.old().pos(exp_i) = conf.old().pos(i) + shift;
	}
      }
    }
  }

  exp_conf.current().box(0) = sim.param().multicell.x * conf.current().box(0);
  exp_conf.current().box(1) = sim.param().multicell.y * conf.current().box(1);
  exp_conf.current().box(2) = sim.param().multicell.z * conf.current().box(2);
  
  exp_conf.current().energies.zero();
  exp_conf.current().perturbed_energy_derivatives.zero();
  exp_conf.boundary_type = conf.boundary_type;

}


/**
 * calculate nonbonded forces and energies.
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
  
  Pairlist temp;
  temp.resize(topo.num_atoms());
  
  for(unsigned int atom_i = 0; atom_i < topo.num_atoms(); ++atom_i){
    
    for(int i=0; i < m_set_size; ++i){

      assert (m_nonbonded_set.size() > i);
      assert (m_nonbonded_set[i]->pairlist().size() > atom_i);
      
      for(unsigned int atom_j = 0;
	  atom_j < m_nonbonded_set[i]->pairlist()[atom_i].size();
	  ++atom_j){
	
	assert(temp.size() > atom_i);
	assert(temp.size() > m_nonbonded_set[i]->pairlist()[atom_i][atom_j]);

	if (m_nonbonded_set[i]->pairlist()[atom_i][atom_j] < atom_i)
	  temp[m_nonbonded_set[i]->pairlist()[atom_i][atom_j]].push_back(atom_i);
	else
	  temp[atom_i].push_back(m_nonbonded_set[i]->pairlist()[atom_i][atom_j]);
      }
      
    }
  }
  
  os << temp << std::endl;
  return 0;
}
