/**
 * @file nonbonded_virial_interaction.tcc
 * template methods of Nonbonded_Virial_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::Nonbonded_Virial_Interaction(virial_enum virial)
  : Nonbonded_Interaction<t_simulation, t_pairlist>::Nonbonded_Interaction(),
    m_virial_type(virial)
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::~Nonbonded_Virial_Interaction()
{
  DEBUG(4, "Nonbonded_Virial_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::calculate_interactions(t_simulation &sim)
{
  // prepare for the virial calculation
  m_com_pos.resize(sim.system().pos().size());
  
  simulation::Molecule_Iterator m_it = sim.topology().molecule_begin(),
    m_to = sim.topology().molecule_end();
  
  math::Vec com_pos, com_ekin;
  m_tot_ekin = 0.0;
  
  for( ; m_it != m_to; ++m_it){
    m_it.com(sim.system(), sim.topology().mass(), com_pos, com_ekin);
    m_tot_ekin += com_ekin;
    
    simulation::Atom_Iterator a_it = m_it.begin(),
      a_to = m_it.end();
    
    for( ; a_it != a_to; ++a_it){
      m_com_pos(*a_it) = com_pos;
    }

  }
  
  // do the work...
  Nonbonded_Interaction<t_simulation, t_pairlist>::calculate_interactions(sim);

}

/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::do_interactions(t_simulation &sim, typename t_pairlist::iterator it, 
		  typename t_pairlist::iterator to,
		  typename Nonbonded_Interaction<t_simulation, t_pairlist>::nonbonded_type_enum range)
{
  
  Nonbonded_Interaction<t_simulation, t_pairlist>::do_interactions(sim, it, to, range);
  
}

/**
 * helper function to calculate 1,4 forces and energies.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::do_14_interactions(t_simulation &sim)
{
  
  Nonbonded_Interaction<t_simulation, t_pairlist>::do_14_interactions(sim);
  
}


