/**
 * @file nonbonded_interaction.tcc
 * template methods of Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::Nonbonded_Interaction(t_simulation &sim)
  : Interaction<t_simulation, t_interaction_spec>("NonBonded"),
    Nonbonded_Base(),
    Storage(),
    t_interaction_spec::nonbonded_innerloop_type(*dynamic_cast<Nonbonded_Base *>(this)),
    m_pairlist_algorithm()
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline 
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::~Nonbonded_Interaction()
{
  DEBUG(4, "Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void 
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  // initialize the constants
  if (!sim.steps())
    initialize(sim);

  // need to update pairlist?
  if(!(sim.steps() % sim.nonbonded().update())){
    // create a pairlist
    // zero the longrange forces, energies, virial
    force() = 0.0;
    energies().zero();
    DEBUG(15, "zero the longrange lambda energies");
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
  if (t_interaction_spec::do_virial){
    DEBUG(7, "\tadd long range virial");
    for(size_t i=0; i<3; ++i)
      for(size_t j=0; j<3; ++j)
	sim.system().virial()(i,j) =
	  sim.system().virial()(i,j) + virial()(i,j);
  }
  
  // add 1,4 - interactions
  DEBUG(7, "\t1,4 - interactions");
  do_14_interactions(sim);

  // possibly do the RF contributions due to excluded atoms
  if(sim.nonbonded().RF_exclusion()){
    DEBUG(7, "\tRF excluded interactions and self term");
    do_RF_excluded_interactions(sim);
  }
}

/**
 * add a shortrange interaction
 */
template<typename t_simulation, typename t_interaction_spec>
inline void
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::add_shortrange_pair(t_simulation const & sim, size_t const i, size_t const j)
{
  assert(pairlist().size() > i);
  pairlist()[i].push_back(j);
}

/**
 * add a longrange interaction
 */
template<typename t_simulation, typename t_interaction_spec>
inline void
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::add_longrange_pair(t_simulation & sim, size_t const i, size_t const j)
{
  interaction_innerloop(sim, i, j, *this);
}


//==================================================
// ACCESSORS
//==================================================

/**
 * pairlist accessor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Pairlist & 
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::pairlist()
{
  return m_pairlist;
}
/**
 * const pairlist accessor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Pairlist const &
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::pairlist()const
{
  return m_pairlist;
}
/**
 * perturbed pairlist accessor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Pairlist & 
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::perturbed_pairlist()
{
  return m_perturbed_pairlist;
}
/**
 * const perturbed pairlist accessor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Pairlist const &
interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::perturbed_pairlist()const
{
  return m_perturbed_pairlist;
}

//==================================================
// interaction loops
//==================================================

/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::do_interactions(t_simulation &sim, Pairlist::iterator it, 
		  Pairlist::iterator to)
{  
  DEBUG(7, "\tcalculate interactions");  

  for( ; it != to; ++it){    

    // shortrange, therefore store in simulation.system()
    interaction_innerloop(sim, it.i(), *it, sim.system());

  }
  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::do_14_interactions(t_simulation &sim)
{
  DEBUG(7, "\tcalculate 1,4-interactions");

  std::set<int>::const_iterator it, to;
  
  for(size_t i=0; i<sim.topology().num_solute_atoms(); ++i){
    it = sim.topology().one_four_pair(i).begin();
    to = sim.topology().one_four_pair(i).end();
    
    for( ; it != to; ++it){

      one_four_interaction_innerloop(sim, i, *it);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::do_RF_excluded_interactions(t_simulation &sim)
{
  
  DEBUG(7, "\tcalculate RF excluded interactions");
  
  for(size_t i=0; i<sim.topology().num_solute_atoms(); ++i){
    
    RF_excluded_interaction_innerloop(sim, i);
    
  } // loop over solute atoms

  // Solvent
  simulation::chargegroup_iterator cg_it = sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();
  cg_it += sim.topology().num_solute_chargegroups();
  
  for( ; cg_it != cg_to; ++cg_it){

    RF_solvent_interaction_innerloop(sim, cg_it);

  } // loop over solvent charge groups
}  

/**
 * initialize the arrays
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Nonbonded_Interaction<t_simulation, t_interaction_spec>
::initialize(t_simulation &sim)
{
  DEBUG(15, "nonbonded_interaction::initialize");
  
  Nonbonded_Base::initialize(sim);

  force().resize(sim.system().force().size());

  energies().resize(sim.system().energies().bond_energy.size(),
		    sim.system().energies().kinetic_energy.size());

  lambda_energies().resize(sim.system().lambda_energies().bond_energy.size(),
			   sim.system().lambda_energies().kinetic_energy.size());
  
}
