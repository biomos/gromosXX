/**
 * @file topology.cc
 * methods definition
 */

#include <stdheader.h>

double topology_ver = 0.30;

namespace simulation
{
  class Multibath;
}


#include <util/debug.h>
#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/simulation.h>

#include <interaction/interaction_types.h>

#undef MODULE
#undef SUBMODULE
#define MODULE topology
#define SUBMODULE topology

/**
 * Constructor
 */
topology::Topology::Topology()
  : m_mass(0),
    m_charge(0),
    m_num_solute_chargegroups(0),
    m_num_solute_molecules(0),
    m_multicell_topo(NULL),
    m_polarizability(0),
    m_coscharge(0),
    m_dadl(0),
    m_damping_level(0),
    m_damping_power(0)
{
  m_chargegroup.push_back(0);
  m_molecule.push_back(0);
}

topology::Topology::~Topology()
{
  delete m_multicell_topo;
}

/**
 * copy (multiply) constructor
 */
topology::Topology::Topology(topology::Topology const & topo, int mul_solute, int mul_solvent)
  : m_multicell_topo(NULL)
{

  DEBUG(7, "multiplying topology");
  
  if (mul_solvent == -1) mul_solvent = mul_solute;
  
  const int num_solute = topo.num_solute_atoms();
  //assert(num_solute >= 0 && topo.m_is_perturbed.size() == unsigned(num_solute));
  
  m_is_perturbed.clear();
  m_is_polarizable.clear();
  m_iac.clear();
  m_mass.clear();
  m_charge.clear();
  m_polarizability.clear();
  m_coscharge.clear();
  m_dadl.clear();
  m_damping_level.clear();
  m_damping_power.clear();
  m_exclusion.clear();
  m_one_four_pair.clear();
  m_all_exclusion.clear();
  m_atom_energy_group.clear();
  m_chargegroup.clear();
  m_chargegroup.push_back(0);
  
  DEBUG(10, "solute chargegrous = " << topo.num_solute_chargegroups());
  
  m_num_solute_chargegroups = topo.num_solute_chargegroups() * mul_solute;
  m_num_solute_molecules = m_num_solute_molecules * mul_solute;

  m_molecule.clear();
  m_molecule.push_back(0);
  
  solute().bonds().clear();
  solute().angles().clear();
  solute().improper_dihedrals().clear();
  solute().dihedrals().clear();

  DEBUG(8, "\tmultiplying solute");
  
  for(int m=0; m<mul_solute; ++m){
    DEBUG(10, "\tmul " << m);

    m_is_perturbed.insert(m_is_perturbed.end(), topo.m_is_perturbed.begin(), topo.m_is_perturbed.end());
    m_is_polarizable.insert(m_is_polarizable.end(), topo.m_is_polarizable.begin(), topo.m_is_polarizable.end());

    for(int i=0; i<num_solute; ++i){
      m_iac.push_back(topo.m_iac[i]);
      m_mass.push_back(topo.m_mass[i]);
      m_charge.push_back(topo.m_charge[i]);
      m_polarizability.push_back(topo.m_polarizability[i]);
      m_coscharge.push_back(topo.m_coscharge[i]);
      m_dadl.push_back(topo.m_dadl[i]);
      m_damping_level.push_back(topo.m_damping_level[i]);
      m_damping_power.push_back(topo.m_damping_power[i]);
      m_exclusion.push_back(topo.m_exclusion[i]);
      m_one_four_pair.push_back(topo.m_one_four_pair[i]);
      m_all_exclusion.push_back(topo.m_all_exclusion[i]);
      m_atom_energy_group.push_back(topo.m_atom_energy_group[i]);
      
      solute().add_atom(topo.solute().atom(i).name, topo.solute().atom(i).residue_nr);

    }

    DEBUG(10, "\tcg");
    for(unsigned int i=1; i<=topo.num_solute_chargegroups(); ++i){
      m_chargegroup.push_back(topo.m_chargegroup[i] + m * num_solute);
    }

    DEBUG(10, "\tmol");
    for(unsigned int i=1; i<topo.molecules().size() - topo.num_solvent_molecules(0); ++i){
      m_molecule.push_back(topo.molecules()[i] + m * num_solute);
    }

    DEBUG(10, "\tbonds");
    for(unsigned int i=0; i<topo.solute().bonds().size(); ++i){
      solute().bonds().push_back
	(two_body_term_struct(topo.solute().bonds()[i].i + m * num_solute,
			      topo.solute().bonds()[i].j + m * num_solute,
			      topo.solute().bonds()[i].type));
    }

    DEBUG(10, "\tdistance constraints");
    for(unsigned int i=0; i<topo.solute().distance_constraints().size(); ++i){
      solute().add_distance_constraint
	(two_body_term_struct(topo.solute().distance_constraints()[i].i + m * num_solute,
			      topo.solute().distance_constraints()[i].j + m * num_solute,
			      topo.solute().distance_constraints()[i].type));
    }

    DEBUG(10, "\tangles");
    for(unsigned int i=0; i<topo.solute().angles().size(); ++i){
      solute().angles().push_back
	(three_body_term_struct(topo.solute().angles()[i].i + m * num_solute,
				topo.solute().angles()[i].j + m * num_solute,
				topo.solute().angles()[i].k + m * num_solute,
				topo.solute().angles()[i].type));
    }

    DEBUG(10, "\timps");
    for(unsigned int i=0; i<topo.solute().improper_dihedrals().size(); ++i){
      solute().improper_dihedrals().push_back
	(four_body_term_struct(topo.solute().improper_dihedrals()[i].i + m * num_solute,
			       topo.solute().improper_dihedrals()[i].j + m * num_solute,
			       topo.solute().improper_dihedrals()[i].k + m * num_solute,
			       topo.solute().improper_dihedrals()[i].l + m * num_solute,
			       topo.solute().improper_dihedrals()[i].type));
    }

    DEBUG(10, "\tdihedrals");
    for(unsigned int i=0; i<topo.solute().dihedrals().size(); ++i){
      solute().dihedrals().push_back
	(four_body_term_struct(topo.solute().dihedrals()[i].i + m * num_solute,
			       topo.solute().dihedrals()[i].j + m * num_solute,
			       topo.solute().dihedrals()[i].k + m * num_solute,
			       topo.solute().dihedrals()[i].l + m * num_solute,
			       topo.solute().dihedrals()[i].type));
    }

    // perturbed solute
    DEBUG(10, "\tperturbed solute: MISSING");
    
  }
  
  DEBUG(10, "\tegroups");
  m_energy_group = topo.m_energy_group;

  DEBUG(10, "\tperturbation param");
  m_lambda = topo.m_lambda;
  m_old_lambda = topo.m_old_lambda;
  m_lambda_exp = topo.m_lambda_exp;
  m_energy_group_scaling = topo.m_energy_group_scaling;
  m_energy_group_lambdadep = topo.m_energy_group_lambdadep;
  m_lambda_prime = topo.m_lambda_prime;
  m_lambda_prime_derivative = topo.m_lambda_prime_derivative;
  m_perturbed_energy_derivative_alpha =   topo.m_perturbed_energy_derivative_alpha;
  
  DEBUG(10, "\tspecial");
  m_position_restraint = topo.m_position_restraint;
  m_jvalue_restraint = topo.m_jvalue_restraint;

  DEBUG(8, "\tnum_solute_atoms()=" << num_solute_atoms());

  if (topo.num_solvent_atoms()) {
    // solvent
    DEBUG(8, "\tmultiplying solvent");

    assert(topo.num_solvents() == 1);
    add_solvent(topo.solvent(0));
    solvate(0, topo.num_solvent_atoms() * mul_solvent / topo.solvent(0).num_atoms());
    
    for(int m = 0; m < mul_solvent; ++m)
      for(unsigned int s = topo.num_solute_atoms(); s < topo.num_atoms(); ++s){
        m_atom_energy_group.push_back(topo.m_atom_energy_group[s]);
    }
  
    DEBUG(8, "\tnum_atoms()=" << num_atoms());
  }
}

/**
 * set the capacity of solute atoms by resizeing
 * the apropriate arrays.
 */
void topology::Topology::resize(unsigned int const atoms)
{
  DEBUG(8, "resizing topology for " << atoms << " atoms!");
  
  // if you want to shrink, first change num_atoms...
  assert(atoms >= num_solute_atoms() + num_solvent_atoms());

  // m_mass.resizeAndPreserve(atoms);
  m_mass.resize(atoms);
  // m_charge.resizeAndPreserve(atoms);
  m_charge.resize(atoms);
  m_polarizability.resize(atoms);
  m_coscharge.resize(atoms);
  m_dadl.resize(atoms);
  m_damping_level.resize(atoms);
  m_damping_power.resize(atoms);
  m_exclusion.resize(atoms);
  m_one_four_pair.resize(atoms);
  m_all_exclusion.resize(atoms);
  
  m_iac.resize(atoms);
  // chargegroups???
  m_is_perturbed.resize(atoms, false);
  m_is_polarizable.resize(atoms, false);
}

void topology::Topology::init(simulation::Simulation const & sim, std::ostream & os, bool quiet)
{
  // std::cerr << "topology::init" << std::endl;
  
  if (!m_molecule.size()){
    m_molecule.push_back(0);
    // do it gromos++ like and determine by bonds???
    m_molecule.push_back(num_solute_atoms());
  }

  if (!m_energy_group.size()){
    m_energy_group.push_back(num_solute_atoms()-1);
    for(unsigned int i=0; i<num_solute_atoms(); ++i)
      m_atom_energy_group.push_back(0);
  }

  if (sim.param().multicell.multicell){
    
    const int num = sim.param().multicell.x * sim.param().multicell.y * sim.param().multicell.z;
    simulation::Simulation multi_sim = sim;
    multi_sim.param().multicell.multicell = false;
    m_multicell_topo = new Topology(*this, num);
    m_multicell_topo->init(multi_sim, os, true);
    m_multicell_topo->check_state();

    if (!quiet){
      os << "\tMULTICELL\n" << "\t\t"
		<< std::setw(10) << sim.param().multicell.x
		<< std::setw(10) << sim.param().multicell.y
		<< std::setw(10) << sim.param().multicell.z
		<< "\n\t\t"
		<< "molecules : " << m_multicell_topo->molecules().size() - 1
		<< "\n\t\t"
		<< "atoms     : " << m_multicell_topo->num_atoms()
		<< "\n\tEND\n";
    }
  }

  // add chargegroup exclusions (a clever technique to improve pairlisting...)
  m_chargegroup_exclusion.resize(num_solute_chargegroups());

  for(size_t cg1=0; cg1<num_solute_chargegroups(); ++cg1){

    m_chargegroup_exclusion[cg1].insert(cg1);
    
    for(size_t cg2=cg1+1; cg2 < num_solute_chargegroups(); ++cg2){

      // std::cerr << "\tchecking cg1=" << cg1 << " cg2=" << cg2 << std::endl;
      
      for(int at1 = m_chargegroup[cg1]; at1 < m_chargegroup[cg1+1]; ++at1){

	std::set<int>::iterator ex = m_all_exclusion[at1].begin(),
	  ex_to = m_all_exclusion[at1].end();
	for( ; ex != ex_to; ++ex){

	  // std::cerr << "cg2: " << m_chargegroup[cg2] << " - " << m_chargegroup[cg2+1]
	  // << " ex=" << *ex << std::endl;
	  
	  if (m_chargegroup[cg2] <= *ex &&
	      m_chargegroup[cg2+1] > *ex){
	  
	    // std::cerr << "cg1=" << cg1 << " cg2=" << cg2 << " : excluded!" << std::endl;

	    m_chargegroup_exclusion[cg1].insert(cg2);
	    m_chargegroup_exclusion[cg2].insert(cg1);

	  } // exclude cg
	} // exclusions of at1
      } // at1 in cg1
    } // cg2
  } // cg1

  if (!sim.param().force.nonbonded_crf){
    for(unsigned int i=0; i<m_charge.size(); ++i){
      m_charge(i) = 0.0;
    }
  }

  if (!sim.param().force.nonbonded_vdw){
    if (m_atom_name.find("DUM") == m_atom_name.end()){
      io::messages.add("no dummy atomtype (DUM) found in topology",
		       "topology",
		       io::message::error);
    }
    else{
      int dum = m_atom_name["DUM"];
      for(unsigned int i=0; i<m_iac.size(); ++i)
	m_iac[i] = dum;
    }
  }

}

/**
 * add a solute atom to the topology.
 * if the arrays are too small they will be increased.
 * if adding multiple solute atoms, first call solute_atoms_capacity...
 */
void topology::Topology
::add_solute_atom(std::string name, int residue_nr,
		  int iac, double mass,
		  double charge, bool chargegroup,
		  std::set<int> exclusions,
		  std::set<int> one_four_pairs)
{

  if (unsigned(m_mass.size()) < num_solute_atoms() + 1){
    resize(num_solute_atoms()+1);
  }
  
  Topology::mass()(num_solute_atoms()) = mass;
  Topology::charge()(num_solute_atoms()) = charge;
  
  Topology::polarizability()(num_solute_atoms()) = 0.0;
  Topology::coscharge()(num_solute_atoms()) = 0.0;
  Topology::dadl()(num_solute_atoms()) = 0.0;
  Topology::damping_level()(num_solute_atoms()) = 0.0;
  Topology::damping_power()(num_solute_atoms()) = 0.0;

  if (chargegroup){
    m_chargegroup.push_back(num_solute_atoms()+1);
    ++m_num_solute_chargegroups;
  }
  
  DEBUG(15, "iac[" << num_solute_atoms() << "] = " << iac);

  m_iac[num_solute_atoms()] = iac;

  m_exclusion[num_solute_atoms()] = exclusions;
  m_one_four_pair[num_solute_atoms()] = one_four_pairs;
  
  std::set_union(exclusions.begin(), exclusions.end(),
		 one_four_pairs.begin(), one_four_pairs.end(),
		 std::inserter(m_all_exclusion[num_solute_atoms()], 
			       m_all_exclusion[num_solute_atoms()].end())
		 );

  // this increases num_solute_atoms()
  solute().add_atom(name, residue_nr);
}

/**
 * add solvent molecules to the simulation (system).
 */
void topology::Topology::solvate(unsigned int solv, unsigned int num_molecules)
{
  // only add in the correct order!
  assert(solv == m_num_solvent_atoms.size());
  assert(solv < m_solvent.size());

  int n = num_solute_atoms() + num_solvent_atoms();

  m_num_solvent_molecules.push_back(num_molecules);
  m_num_solvent_atoms.push_back(num_molecules * m_solvent[solv].num_atoms());
  
  DEBUG(5, "solvate: solvent atoms: " << num_solvent_atoms());
  DEBUG(10, "solvate: total atoms: " << num_solute_atoms() + num_solvent_atoms());
  
  resize(num_solute_atoms() + num_solvent_atoms());

  // add to iac, mass, charge
  for(unsigned int i=0; i<num_molecules; ++i){
    for(unsigned int j=0; j<m_solvent[solv].num_atoms(); ++j, ++n){

      DEBUG(15, "iac[" << n << "]=" << m_solvent[solv].atom(j).iac);
      DEBUG(15, "charge[" << n << "]=" << m_solvent[solv].atom(j).charge);
      
      m_iac[n] = m_solvent[solv].atom(j).iac;
      m_mass(n) = m_solvent[solv].atom(j).mass;
      m_charge(n) = m_solvent[solv].atom(j).charge;
      
      m_polarizability[n] = m_solvent[solv].atom(j).polarizability;
      m_coscharge(n) = m_solvent[solv].atom(j).coscharge;
      m_damping_level(n) = m_solvent[solv].atom(j).damping_level;
      m_damping_power(n) = m_solvent[solv].atom(j).damping_power;
      // no exclusions or 1-4 interactions for solvent ?!

    }

    // add to the chargegroups
    DEBUG(8, "solvent cg: " << n);
    m_chargegroup.push_back(n);

    // and to the molecules
    m_molecule.push_back(n);

  }
    
}

void topology::Topology::resolvate(unsigned int solv, unsigned int num_molecules)
{
  if (solv != 0){
    io::messages.add("resolvation for multiple solvents not implemented",
		     "topology",
		     io::message::error);
    return;
  }
  
  assert(m_num_solvent_atoms.size() == 1);
  assert(m_solvent.size() == 1);

  int n = num_solute_atoms();

  m_num_solvent_molecules[m_num_solvent_molecules.size() - 1] = num_molecules;
  m_num_solvent_atoms[m_num_solvent_atoms.size() - 1] = num_molecules * m_solvent[solv].num_atoms();
  
  DEBUG(5, "solvate: solvent atoms: " << num_solvent_atoms());
  DEBUG(10, "solvate: total atoms: " << num_solute_atoms() + num_solvent_atoms());
  
  DEBUG(8, "resizing for solute: " << num_solute_atoms() << " and solvent: "
	<< num_solvent_atoms() << " (total: " << num_solute_atoms() + num_solvent_atoms());
  resize(num_solute_atoms() + num_solvent_atoms());

  m_chargegroup.erase(m_chargegroup.begin() + m_num_solute_chargegroups + 1, m_chargegroup.end());
  m_molecule.erase(m_molecule.begin() + m_num_solute_molecules + 1, m_molecule.end());

  DEBUG(9, "solute chargegroups: " << m_num_solute_chargegroups
	<< " size: " << m_chargegroup.size());
  DEBUG(9, "solute molecules: " << m_num_solute_molecules
	<< " size: " << m_molecule.size());

  // add to iac, mass, charge
  for(unsigned int i=0; i<num_molecules; ++i){
    for(unsigned int j=0; j<m_solvent[solv].num_atoms(); ++j, ++n){

      DEBUG(15, "iac[" << n << "]=" << m_solvent[solv].atom(j).iac);
      DEBUG(15, "charge[" << n << "]=" << m_solvent[solv].atom(j).charge);
      
      m_iac[n] = m_solvent[solv].atom(j).iac;
      m_mass(n) = m_solvent[solv].atom(j).mass;
      m_charge(n) = m_solvent[solv].atom(j).charge;
      
      m_polarizability(n) = m_solvent[solv].atom(j).polarizability;
      m_coscharge(n) = m_solvent[solv].atom(j).coscharge;
      m_damping_level(n) = m_solvent[solv].atom(j).damping_level;
      m_damping_power(n) = m_solvent[solv].atom(j).damping_power;
      
      // no exclusions or 1-4 interactions for solvent
    }

    // add to the chargegroups
    DEBUG(8, "solvent cg: " << n);
    m_chargegroup.push_back(n);

    // and to the molecules
    DEBUG(8, "solvent mol: " << n);
    m_molecule.push_back(n);
  }
}

/**
 * total number of solvent atoms.
 */
unsigned int topology::Topology::num_solvent_atoms()const
{
  unsigned int n = 0;
  for(std::vector<unsigned int>::const_iterator it = m_num_solvent_atoms.begin(),
	to = m_num_solvent_atoms.end();
      it != to; ++it)
    n += *it;
  return n;
}

/**
 * calculate constraint degrees of freedom.
 */
void
topology::Topology::
calculate_constraint_dof(simulation::Multibath &multibath,
			 bool rottrans_constraints)const
{
  // substract constraints
  {
    std::vector<two_body_term_struct>::const_iterator 
      c_it = solute().distance_constraints().begin(),
      c_to = solute().distance_constraints().end();
    
    unsigned int com_bath_i, ir_bath_i, com_bath_j, ir_bath_j;
    
    for( ; c_it != c_to; ++c_it){
      
      DEBUG(10, "Constraint: " << c_it->i << " - " << c_it->j);
      multibath.in_bath(c_it->i, com_bath_i, ir_bath_i);
      multibath.in_bath(c_it->j, com_bath_j, ir_bath_j);
      
      multibath[ir_bath_i].dof -= 0.5;
      multibath[ir_bath_j].dof -= 0.5;
      
      multibath[ir_bath_i].ir_dof -= 0.5;
      multibath[ir_bath_j].ir_dof -= 0.5;
      
      multibath[ir_bath_i].solute_constr_dof += 0.5;
      multibath[ir_bath_j].solute_constr_dof += 0.5;
      
    }
    
    for(unsigned int i=0; i<multibath.size(); ++i){
      DEBUG(7, "dof           " << multibath[i].dof);
      DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
    }
    
    // solvent constraints
    int index = num_solute_atoms();
    for(unsigned int s=0; s < num_solvents(); ++s){
      
      for(unsigned int m=0; m < num_solvent_molecules(s); ++m){
	
	c_it = solvent(s).distance_constraints().begin();
	c_to = solvent(s).distance_constraints().end();
	
	for( ; c_it != c_to; ++c_it){
	  
	  multibath.in_bath(c_it->i + index, com_bath_i, ir_bath_i);
	  multibath.in_bath(c_it->j + index, com_bath_j, ir_bath_j);
	  
	  multibath[ir_bath_i].dof -= 0.5;
	  multibath[ir_bath_j].dof -= 0.5;
	  
	  multibath[ir_bath_i].ir_dof -= 0.5;
	  multibath[ir_bath_j].ir_dof -= 0.5;
	  
	  multibath[ir_bath_i].solvent_constr_dof += 0.5;
	  multibath[ir_bath_j].solvent_constr_dof += 0.5;
	  
	}
	
	index += solvent(s).num_atoms();
	
      }
    }
  }
  
  DEBUG(7, "and the perturbd distance constraints (DOF calc)");
  
  {
    // substract perturbed constraints
    std::vector<perturbed_two_body_term_struct>
      ::const_iterator 
      c_it = perturbed_solute().distance_constraints().begin(),
      c_to = perturbed_solute().distance_constraints().end();
    
    unsigned int com_bath_i, ir_bath_i, com_bath_j, ir_bath_j;
    
    for( ; c_it != c_to; ++c_it){
      
      DEBUG(10, "Constraint: " << c_it->i << " - " << c_it->j);
      multibath.in_bath(c_it->i, com_bath_i, ir_bath_i);
      multibath.in_bath(c_it->j, com_bath_j, ir_bath_j);
      
      multibath[ir_bath_i].dof -= 0.5;
      multibath[ir_bath_j].dof -= 0.5;
      
      multibath[ir_bath_i].ir_dof -= 0.5;
      multibath[ir_bath_j].ir_dof -= 0.5;
      
      multibath[ir_bath_i].solute_constr_dof += 0.5;
      multibath[ir_bath_j].solute_constr_dof += 0.5;
      
    }
  
    for(unsigned int i=0; i<multibath.size(); ++i){
      DEBUG(7, "dof           " << multibath[i].dof);
      DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
    }
  }
  
  if (rottrans_constraints){
    // check whether all solute is in one bath
    if (num_solute_atoms()){
     
      unsigned int ir_bath, ir_bath_0, com_bath, com_bath_0;
      bool ok = true;

      multibath.in_bath(0, com_bath_0, ir_bath_0);
      for(unsigned int i=1; i<num_solute_atoms(); ++i){
	multibath.in_bath(i, com_bath, ir_bath);
	if (com_bath != com_bath_0 || ir_bath != ir_bath_0){
	  io::messages.add("roto-translational constraints: all solute has to be coupled "
			   "to one ir and one com bath", "calc_dof", 
			   io::message::error);
	  ok = false;
	  break;
	}
      }
      if (ok){
	std::cout << "roto-translational constraints: removing 3 dof from com bath "
		  << com_bath_0+1 << " and from ir bath " << ir_bath_0+1 << "\n";
	multibath[com_bath_0].dof -= 3.0;
	multibath[com_bath_0].com_dof -= 3.0;
	multibath[ir_bath_0].dof -= 3.0;
	multibath[ir_bath_0].ir_dof -= 3.0;
      }
    }
  }

  // check whether we have dof in every (coupled) bath
  for(unsigned int i = 0; i < multibath.size(); ++i){
    if (multibath[i].dof <= 0 && multibath[i].tau != -1){
      io::messages.add("removing coupling of bath with 0 dof",
		       "calc_dof",
		       io::message::notice);
      multibath[i].tau = -1;
    }
  }

  DEBUG(10, "end dof calc");

}

void
topology::Topology::update_for_lambda()
{
  for(std::map<unsigned int, topology::Perturbed_Atom>::const_iterator
        it = perturbed_solute().atoms().begin(),
        to = perturbed_solute().atoms().end();
      it != to; ++it){
    mass()(it->second.sequence_number()) = 
      (1-lambda()) * it->second.A_mass() + 
      lambda() * it->second.B_mass();

    DEBUG(8, "mass A : " << it->second.A_mass() << " B : "
          << it->second.B_mass());
    DEBUG(8, "mass(" << it->second.sequence_number()
          << ") = " << mass()(it->second.sequence_number()));
  }

  // this is nowadays done directly in the shake routines
  // perturbed_solute().set_distance_constraints(lambda());
  
}

/**
 * check state
 */
int 
topology::Topology::check_state()const
{
  int result = 0;

  // check that we have an energy group for every atom
  if (m_atom_energy_group.size() != num_atoms()){
    io::messages.add("not every atom has an energy group index",
		     "Topology::check_state", io::message::error);
    ++result;
  }
  for(std::vector<unsigned int>::const_iterator it = m_atom_energy_group.begin(),
	to = m_atom_energy_group.end(); it != to; ++it){
    if (*it >= m_energy_group.size()){
      io::messages.add("energy group index of atom too large",
		       "Topology::check_state", io::message::error);
      ++result;
    }
  }

  return result;
}

namespace topology
{
  /**
   * output information about the topology.
   */
  std::ostream & operator<<(std::ostream &os, Topology &topo)
  {
    os << "a topology";
    return os;
  }
}

