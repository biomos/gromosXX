/**
 * @file nonbonded_interaction.tcc
 * template methods of Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::Nonbonded_Interaction(t_simulation &sim)
  : Interaction<t_simulation>("NonBonded"),
    Nonbonded_Base(),
    t_innerloop(*this, sim.system()),
    m_pairlist(*this)
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::~Nonbonded_Interaction()
{
  DEBUG(4, "Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  // initialize the constants
  if (!sim.steps())
    initialize(sim);

  if(sim.pressure_calculation())
    sim.calculate_mol_com();

  // need to update pairlist?
  DEBUG(7, "steps " << sim.steps() << " upd " << sim.nonbonded().update());

  if(!(sim.steps() % sim.nonbonded().update())){
    // create a pairlist
    DEBUG(7, "\tupdate the parlist");
    m_pairlist.update(sim);
    DEBUG(7, "\tafter update : " << m_pairlist.size());
   
  }

  // calculate forces / energies
  DEBUG(7, "\tshort range");
  // do_interactions(sim, m_pairlist.short_range().begin(),
  do_interactions(sim, m_pairlist.begin(),
		  // m_pairlist.short_range().end(),
		  m_pairlist.end() );
  //		  shortrange);

  
  // add long-range force
  sim.system().force() += m_pairlist.filter().force();
  
  // and long-range energies
  for(size_t i = 0; i < m_pairlist.filter().energies().lj_energy.size(); ++i){
    for(size_t j = 0; j < m_pairlist.filter().energies().lj_energy.size(); ++j){
      sim.system().energies().lj_energy[i][j] += 
	m_pairlist.filter().energies().lj_energy[i][j];
      sim.system().energies().crf_energy[i][j] += 
	m_pairlist.filter().energies().crf_energy[i][j];
    }
  }

  // add longrange virial
  if (sim.pressure_calculation()){
    for(size_t i=0; i<3; ++i)
      for(size_t j=0; j<3; ++j)
	sim.system().virial()(i,j) =
	  -0.5 * (sim.system().virial()(i,j) + m_pairlist.filter().virial()(i,j));
  }
  
  // add 1,4 - interactions
  do_14_interactions(sim);

  // possibly do the RF contributions due to excluded atoms
  if(sim.nonbonded().RF_exclusion())
    do_RF_excluded_interactions(sim);

}

/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::do_interactions(t_simulation &sim, typename t_pairlist::iterator it, 
		  typename t_pairlist::iterator to)
{  
  DEBUG(7, "\tcalculate interactions");  

  for( ; it != to; ++it){    

    interaction_inner_loop(sim, it.i(), *it);

  }
  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::do_14_interactions(t_simulation &sim)
{
  DEBUG(7, "\tcalculate 1,4-interactions");

  std::set<int>::const_iterator it, to;
  
  for(size_t i=0; i<sim.topology().num_solute_atoms(); ++i){
    it = sim.topology().one_four_pair(i).begin();
    to = sim.topology().one_four_pair(i).end();
    
    for( ; it != to; ++it){

      one_four_interaction_inner_loop(sim, i, *it);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::do_RF_excluded_interactions(t_simulation &sim)
{
  
  DEBUG(7, "\tcalculate RF excluded interactions");
  
  for(size_t i=0; i<sim.topology().num_solute_atoms(); ++i){
    
    RF_excluded_interaction_inner_loop(sim, i);
    
  } // loop over solute atoms

  // Solvent
  simulation::chargegroup_iterator cg_it = sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();
  cg_it += sim.topology().num_solute_chargegroups();
  
  for( ; cg_it != cg_to; ++cg_it){

    RF_solvent_interaction_inner_loop(sim, cg_it);

  } // loop over solvent charge groups
}  

/**
 * pairlist accessor
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
t_pairlist & 
interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::pairlist()
{
  return m_pairlist;
}
