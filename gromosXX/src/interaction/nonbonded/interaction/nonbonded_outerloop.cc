/**
 * @file nonbonded_interaction.cc
 * template methods of Nonbonded_Interaction.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <math/periodicity.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/nonbonded_innerloop.h>

#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>

#include <interaction/nonbonded/interaction_spec.h>

#include <util/debug.h>
#include <util/template_split.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Nonbonded_Outerloop
::Nonbonded_Outerloop(Nonbonded_Parameter &nbp)
  : m_param(nbp)
{
}

//==================================================
// interaction loops
//==================================================

void interaction::Nonbonded_Outerloop
::lj_crf_outerloop(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim, 
		   Pairlist const & pairlist,
		   Storage & storage)
{
  SPLIT_INNERLOOP(_lj_crf_outerloop, topo, conf, sim, pairlist, storage);
}


/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_lj_crf_outerloop(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
		    Pairlist const & pairlist,
		    Storage & storage)
{  
  DEBUG(7, "\tcalculate interactions");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  /*
    variables for a OMP parallelizable loop.
    outer index has to be integer...
  */
  std::vector<unsigned int>::const_iterator j_it, j_to;
  unsigned int i;
  unsigned int size_i = unsigned(pairlist.size());

  //**************************
  // the Bekker implementation
  //**************************
  if (t_interaction_spec::do_bekker){

    periodicity.recalc_shift_vectors();

    int pc;
    unsigned int j;
    // translate the atom j
    DEBUG(9, "nonbonded_interaction: grid based pairlist");

    for(i=0; i < size_i; ++i){

      for(j_it = pairlist[i].begin(),
	    j_to = pairlist[i].end();
	  j_it != j_to;
	  ++j_it){
      
	pc = (*j_it >> 26);
	j = (*j_it & 67108863);
      
	DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << j
	      << " pc " << pc);
      
	innerloop.lj_crf_innerloop(topo, conf, i, j, storage, periodicity, pc);
      }
      
    }

  }
  //*************************
  // standard implementation
  //*************************
  else{ // no grid based pairlist

    DEBUG(9, "nonbonded_interaction: no grid based pairlist");
    for(i=0; i < size_i; ++i){
    
      for(j_it = pairlist[i].begin(),
	    j_to = pairlist[i].end();
	  j_it != j_to;
	  ++j_it){

	DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << *j_it);
	// printf("nb pair %d - %d\n", i, *j_it);
	
	// shortrange, therefore store in simulation.system()
	innerloop.lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);

	// storage.energies.bond_energy[0] += *j_it;

      }
      
    }
  }
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
void interaction::Nonbonded_Outerloop
::one_four_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage)
{
  SPLIT_INNERLOOP(_one_four_outerloop, topo, conf, sim, storage);
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_one_four_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage)
{
  DEBUG(7, "\tcalculate 1,4-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  
  std::set<int>::const_iterator it, to;
  unsigned int const num_solute_atoms = topo.num_solute_atoms();
  unsigned int i;
  
  for(i=0; i < num_solute_atoms; ++i){
    it = topo.one_four_pair(i).begin();
    to = topo.one_four_pair(i).end();
    
    for( ; it != to; ++it){

      innerloop.one_four_interaction_innerloop(topo, conf, i, *it, periodicity);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
void interaction::Nonbonded_Outerloop
::RF_excluded_outerloop(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			Storage & storage)
{
  SPLIT_INNERLOOP(_RF_excluded_outerloop, topo, conf, sim, storage);  
}

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_RF_excluded_outerloop(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			Storage & storage)
{
  
  DEBUG(7, "\tcalculate RF excluded interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  
  int i, num_solute_atoms = topo.num_solute_atoms();
  
  for(i=0; i<num_solute_atoms; ++i){
    
    innerloop.RF_excluded_interaction_innerloop(topo, conf, i, periodicity);
    
  } // loop over solute atoms

  // Solvent
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();
  cg_it += topo.num_solute_chargegroups();

  for(; cg_it != cg_to; ++cg_it){

    innerloop.RF_solvent_interaction_innerloop(topo, conf, cg_it, periodicity);
    ++cg_it;

  } // loop over solvent charge groups
}  

/**
 * calculate the interaction for a given atom pair.
 * SLOW! as it has to create the periodicity...
 */
int interaction::Nonbonded_Outerloop::calculate_interaction
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 unsigned int atom_i, unsigned int atom_j,
 math::Vec & force, 
 double &e_lj, double &e_crf
 )
{
  SPLIT_INNERLOOP(_calculate_interaction, topo, conf, sim, atom_i, atom_j, force, e_lj, e_crf);
  return 0;
}

template<typename t_interaction_spec>
int interaction::Nonbonded_Outerloop
::_calculate_interaction(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 unsigned int atom_i, unsigned int atom_j,
			 math::Vec & force,
			 double & e_lj, double & e_crf)
{
  math::Vec r;
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  
  Nonbonded_Term term;
  term.init(sim);

  const lj_parameter_struct &lj = 
    m_param.lj_parameter(topo.iac(atom_i),
			 topo.iac(atom_j));
  
  periodicity.nearest_image(conf.current().pos(atom_i), conf.current().pos(atom_j), r);

  double f;
  term.lj_crf_interaction(r, lj.c6, lj.c12,
			  topo.charge()(atom_i) * topo.charge()(atom_j),
			  f, e_lj, e_crf);
  force = f * r;

  return 0;
}


/**
 * calculate the hessian for a given atom.
 * this will be VERY SLOW !
 */
int interaction::Nonbonded_Outerloop
::calculate_hessian(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    unsigned int atom_i, unsigned int atom_j,
		    math::Matrix & hessian,
		    Pairlist const & pairlist){

  SPLIT_INNERLOOP(_calculate_hessian, topo, conf, sim, atom_i, atom_j, hessian, pairlist);
  return 0;
}

template<typename t_interaction_spec>
int interaction::Nonbonded_Outerloop
::_calculate_hessian(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     unsigned int atom_i, unsigned int atom_j,
		     math::Matrix & hessian,
		     Pairlist const & pairlist){
  
  hessian = 0.0;
  
  // loop over the pairlist

  //*************************
  // standard implementation
  //*************************

  std::vector<unsigned int>::const_iterator j_it, j_to;

  math::Vec r;
  math::Matrix h;
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  
  Nonbonded_Term term;
  term.init(sim);
  
  assert(pairlist.size() > atom_i);
  
  for(j_it = pairlist[atom_i].begin(),
	j_to = pairlist[atom_i].end();
      j_it != j_to;
      ++j_it){

    DEBUG(8, "\thessian: checking pairlist[" << atom_i << "] : " << *j_it);

    if (*j_it != atom_j) continue;
    DEBUG(8, "\thessian pair in pairlist: " << atom_i << " - " << atom_j);

    periodicity.nearest_image(conf.current().pos(atom_i),
			      conf.current().pos(atom_j),
			      r);
      
    const lj_parameter_struct &lj = 
      m_param.lj_parameter(topo.iac(atom_i),
			   topo.iac(atom_j));
    
    term.lj_crf_hessian(r,
			lj.c6, lj.c12, 
			topo.charge()(atom_i) * topo.charge()(atom_j),
			h);

    for(unsigned int d1=0; d1 < 3; ++d1)
      for(unsigned int d2=0; d2 < 3; ++d2)
	hessian(d1,d2) += h(d1,d2);
  }

  // and the other way round
  assert(pairlist.size() > atom_j);
  for(j_it = pairlist[atom_j].begin(),
	j_to = pairlist[atom_j].end();
      j_it != j_to;
      ++j_it){
    
    if (*j_it != atom_i) continue;
    DEBUG(9, "\thessian pair in pairlist: " << atom_j << " - " << atom_i);

    periodicity.nearest_image(conf.current().pos(atom_i),
			      conf.current().pos(atom_j),
			      r);
      
    const lj_parameter_struct &lj = 
      m_param.lj_parameter(topo.iac(atom_i),
			   topo.iac(atom_j));
    
    term.lj_crf_hessian(r,
			lj.c6, lj.c12, 
			topo.charge()(atom_i) * topo.charge()(atom_j),
			h);

    for(unsigned int d1=0; d1 < 3; ++d1)
      for(unsigned int d2=0; d2 < 3; ++d2)
	hessian(d1,d2) += h(d1,d2);
  }
  
  return 0;
}
