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
template<typename t_simulation, typename t_pairlist>
inline interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::Nonbonded_Interaction(t_simulation &sim)
  : Interaction<t_simulation>("NonBonded"),
    Nonbonded_Base(),
    Nonbonded_Inner_Loop<t_simulation, typename
			 t_simulation::system_type>(*this, sim.system()),
    m_pairlist(*this)
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::~Nonbonded_Interaction()
{
  DEBUG(4, "Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  // initialize the constants
  if (!sim.steps())
    initialize(sim);

  // need to update pairlist?
  DEBUG(7, "steps " << sim.steps() << " upd " << sim.nonbonded().update());

  if(!(sim.steps() % sim.nonbonded().update())){
    // create a pairlist
    DEBUG(7, "\tupdate the parlist");
    m_pairlist.update(sim);
    DEBUG(7, "\tafter update : " << m_pairlist.size());
    
    /*
    // recalc long-range forces
    DEBUG(7, "\tlong range");
    m_longrange_force.resize(sim.system().force().size());
    m_longrange_force = 0.0;

    m_longrange_energy.resize(sim.system().energies().bond_energy.size());
    
    do_interactions(sim, m_pairlist.long_range().begin(),
		    m_pairlist.long_range().end(),
		    longrange);
    */

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
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
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
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
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
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::do_RF_excluded_interactions(t_simulation &sim)
{
  
  math::Vec r, f;
  double e_crf;
  std::cout.precision(10);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  
  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();

  DEBUG(7, "\tcalculate RF excluded interactions");
  
  for(size_t i=0; i<sim.topology().num_solute_atoms(); ++i){

    std::set<int>::const_iterator it, to;
    it = sim.topology().exclusion(i).begin();
    to = sim.topology().exclusion(i).end();
    
    DEBUG(11, "\tself-term " << i );
    r=0;
    
    // this will only contribute in the energy, the force should be zero.
    rf_interaction(r,sim.topology().charge()(i) * sim.topology().charge()(i),
		   f, e_crf);
    sim.system().energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(i)] += 0.5 * e_crf;
    DEBUG(11, "\tcontribution " << 0.5*e_crf);
    
    for( ; it != to; ++it){
      
      DEBUG(11, "\texcluded pair " << i << " - " << *it);
      
      sim.system().periodicity().nearest_image(pos(i), pos(*it), r);
      
      
      rf_interaction(r, sim.topology().charge()(i) * 
		     sim.topology().charge()(*it),
		     f, e_crf);
      
      force(i) += f;
      force(*it) -= f;
      
      // energy
      sim.system().energies().crf_energy[sim.topology().atom_energy_group(i)]
	[sim.topology().atom_energy_group(*it)] += e_crf;
      DEBUG(11, "\tcontribution " << e_crf);
      
    } // loop over excluded pairs
    
    
  } // loop over solute atoms

  // Solvent
  simulation::chargegroup_iterator cg_it = sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();
  cg_it += sim.topology().num_solute_chargegroups();
  
  for( ; cg_it != cg_to; ++cg_it){

    // loop over the atoms
    simulation::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();

    for ( ; at_it != at_to; ++at_it){
      DEBUG(11, "\tsolvent self term " << *at_it);
      // no solvent self term. The distance dependent part and the forces
      // are zero. The distance independent part should add up to zero 
      // for the energies and is left out.

      for(simulation::Atom_Iterator at2_it=at_it+1; at2_it!=at_to; ++at2_it){
	
	DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
	sim.system().periodicity().nearest_image(pos(*at_it), 
						 pos(*at2_it), r);

	// for solvent, we don't calculate internal forces (rigid molecules)
	// and the distance independent parts should go to zero
	e_crf = -sim.topology().charge()(*at_it) * 
	  sim.topology().charge()(*at2_it) * coulomb_constant() * 
	  m_crf_2cut3i * dot(r,r);
	
	// energy
	sim.system().energies().crf_energy
	  [sim.topology().atom_energy_group(*at_it) ]
	  [sim.topology().atom_energy_group(*at2_it)] += e_crf;
      } // loop over at2_it
    } // loop over at_it
  } // loop over solvent charge groups
}  

/**
 * pairlist accessor
 */
template<typename t_simulation, typename t_pairlist>
t_pairlist & 
interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::pairlist()
{
  return m_pairlist;
}
