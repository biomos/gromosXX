/**
 * @file perturbed_nonbonded_interaction.tcc
 * template methods of Perturbed_Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline
interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_nonbonded_spec>
::Perturbed_Nonbonded_Interaction(t_simulation &sim)
  : Nonbonded_Interaction<t_simulation, t_nonbonded_spec>(sim),
    t_nonbonded_spec::perturbation_filter_type(),
    t_nonbonded_spec::perturbed_nonbonded_innerloop_type
      (*dynamic_cast<Nonbonded_Base *>(this))
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline 
interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_nonbonded_spec>
::~Perturbed_Nonbonded_Interaction()
{
  DEBUG(4, "Perturbed_Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline void 
interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_nonbonded_spec>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(4, "(Perturbed) Nonbonded_Interaction::calculate_interactions");

  // initialize the constants
  if (!sim.steps())
    initialize(sim);

  if(sim.pressure_calculation())
    sim.calculate_mol_com();

  // need to update pairlist?
  if(!(sim.steps() % sim.nonbonded().update())){
    // create a pairlist

    // zero the longrange forces, energies, virial
    force() = 0.0;
    energies().zero();
    lambda_energies().zero();
    virial() = 0.0;

    DEBUG(7, "\tupdate the parlist");
    m_pairlist_algorithm.update(sim, *this);
    DEBUG(7, "\tpairlist updated");
  }

  // calculate forces / energies
  DEBUG(7, "\tshort range");
  do_interactions(sim, m_pairlist.begin(),
		  m_pairlist.end() );

  if (t_nonbonded_spec::do_perturbation){
    DEBUG(7, "\tperturbed short range");
    do_perturbed_interactions(sim, m_perturbed_pairlist.begin(),
			      m_perturbed_pairlist.end() );
  }
  
  // add long-range force
  DEBUG(7, "\tadd long range forces and energies");

  sim.system().force() += force();
  
  // and long-range energies
  for(size_t i = 0; i < energies().lj_energy.size(); ++i){
    for(size_t j = 0; j < energies().lj_energy.size(); ++j){
      sim.system().energies().lj_energy[i][j] += 
	energies().lj_energy[i][j];
      sim.system().energies().crf_energy[i][j] += 
	energies().crf_energy[i][j];
    }
  }

  // add longrange virial
  if (t_nonbonded_spec::do_virial){
    DEBUG(7, "\tadd long range virial");
    for(size_t i=0; i<3; ++i)
      for(size_t j=0; j<3; ++j)
	sim.system().virial()(i,j) =
	  -0.5 * (sim.system().virial()(i,j) + virial()(i,j));
  }
  
  // add 1,4 - interactions
  DEBUG(7, "\t1,4 - interactions");
  do_14_interactions(sim);
  if(t_nonbonded_spec::do_perturbation){
    DEBUG(7, "\tperturbed 1,4 - interactions");
    do_perturbed_14_interactions(sim);
  }
  
  // possibly do the RF contributions due to excluded atoms
  if(sim.nonbonded().RF_exclusion()){
    DEBUG(7, "\tRF excluded interactions and self term");
    do_RF_excluded_interactions(sim);
    if(t_nonbonded_spec::do_perturbation){
      DEBUG(7, "\tperturbed RF excluded interactions and self term");
      do_perturbed_RF_excluded_interactions(sim);
    }
  }

  if(t_nonbonded_spec::do_perturbation){
    DEBUG(7, "\tperturbed pairs");
    do_perturbed_pair_interactions(sim);
  }
  
  if (t_nonbonded_spec::do_perturbation){
    // and long-range energy lambda-derivatives
    DEBUG(7, "add long-range lambda-derivatives");

    for(size_t i = 0; 
	i < lambda_energies().lj_energy.size(); ++i){
      for(size_t j = 0; j < lambda_energies().lj_energy.size(); ++j){

	assert(sim.system().lambda_energies().lj_energy.size() > i);
	assert(sim.system().lambda_energies().lj_energy[i].size() > j);
	assert(sim.system().lambda_energies().lj_energy.size() > j);
	assert(sim.system().lambda_energies().lj_energy[j].size() > i);
	
	sim.system().lambda_energies().lj_energy[i][j] += 
	  lambda_energies()
	  .lj_energy[i][j];
	sim.system().lambda_energies().crf_energy[i][j] += 
	  lambda_energies()
	  .crf_energy[i][j];
      }
    }
  } // do perturbed
}

/**
 * add a shortrange interaction
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline void
interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_nonbonded_spec>
::add_shortrange_pair(t_simulation const & sim, size_t const i, size_t const j)
{
  assert(pairlist().size() > i);
  if (perturbed_atom(sim, i))
    perturbed_pairlist()[i].push_back(j);
  else if (perturbed_atom(sim, j))
    perturbed_pairlist()[j].push_back(i);
  else
    pairlist()[i].push_back(j);
}

/**
 * add a longrange interaction
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline void
interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_nonbonded_spec>
::add_longrange_pair(t_simulation & sim, size_t const i, size_t const j)
{
  if (perturbed_atom(sim, i)){
    perturbed_interaction_innerloop(sim, i, j, *this);
  }
  else if (perturbed_atom(sim, j)){
    perturbed_interaction_innerloop(sim, j, i, *this);
  }
  else{
    interaction_innerloop(sim, i, j, *this);
  }
}

//==================================================
// interaction loops
//==================================================

//==================================================
// the perturbed interaction (outer) loops
//==================================================

/**
 * helper function to calculate perturbed forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_nonbonded_spec>
::do_perturbed_interactions(t_simulation &sim,
			    Pairlist::iterator it,
			    Pairlist::iterator to)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  for( ; it != to; ++it){    
    DEBUG(8, "perturbed pair: " << it.i() << " - " << *it);
    perturbed_interaction_innerloop(sim, it.i(), *it, sim.system());
  }

  DEBUG(7, "end of function perturbed nonbonded interaction");  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_nonbonded_spec>
::do_perturbed_14_interactions(t_simulation &sim)
{
  DEBUG(7, "\tcalculate perturbed 1,4-interactions");

  std::set<int>::const_iterator it, to;
  std::map<size_t, simulation::Perturbed_Atom>::const_iterator 
    mit=sim.topology().perturbed_solute().atoms().begin(), 
    mto=sim.topology().perturbed_solute().atoms().end();
  
  for(; mit!=mto; ++mit){
    it = mit->second.one_four_pair().begin();
    to = mit->second.one_four_pair().end();
    
    for( ; it != to; ++it){

      perturbed_one_four_interaction_innerloop
	(sim, mit->second.sequence_number(), *it);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_nonbonded_spec>
::do_perturbed_RF_excluded_interactions(t_simulation &sim)
{

  DEBUG(7, "\tcalculate perturbed excluded RF interactions");

  std::map<size_t, simulation::Perturbed_Atom>::const_iterator
    mit=sim.topology().perturbed_solute().atoms().begin(),
    mto=sim.topology().perturbed_solute().atoms().end();

  DEBUG(9, "\tSize of perturbed atoms " 
	<< sim.topology().perturbed_solute().atoms().size());
  
  for(; mit!=mto; ++mit){
    perturbed_RF_excluded_interaction_innerloop(sim, mit);
  }
}

/**
 * calculate the interactions for the
 * PERTURBED PAIRS
 * (different interaction types in A and in B)
 */
template<typename t_simulation, typename t_nonbonded_spec>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_nonbonded_spec>
::do_perturbed_pair_interactions(t_simulation &sim)
{

  std::vector<simulation::Perturbed_Atompair>::const_iterator
    it = sim.topology().perturbed_solute().atompairs().begin(),
    to = sim.topology().perturbed_solute().atompairs().end();
    
  for(; it != to; ++it){
    perturbed_pair_interaction_innerloop(sim, it);
  }
  
}

