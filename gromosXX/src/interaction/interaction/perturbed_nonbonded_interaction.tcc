/**
 * @file perturbed_nonbonded_interaction.tcc
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
inline interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
						    t_pairlist, t_innerloop>
::Perturbed_Nonbonded_Interaction(t_simulation &sim, interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop> &nonbonded_interaction)
  : Interaction<t_simulation>("Perturbed NonBonded"),
    m_nonbonded_interaction(nonbonded_interaction)
{
  // this should maybe be done somewhere else, but it seems to work
  m_nonbonded_interaction.alpha_lj(sim.topology().alpha_lj());
  m_nonbonded_interaction.alpha_crf(sim.topology().alpha_crf());
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_pairlist, 
						    t_innerloop>
::~Perturbed_Nonbonded_Interaction()
{
  DEBUG(4, "Perturbed_Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
						   t_pairlist, t_innerloop>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(4, "Perturbed_Nonbonded_Interaction::calculate_interactions");

  // calculate forces / energies
  DEBUG(7, "\tshort range");

  do_perturbed_interactions
    (sim, m_nonbonded_interaction.pairlist().perturbed_begin(),
     m_nonbonded_interaction.pairlist().perturbed_end());
  
  // add 1,4 - interactions
  do_perturbed_14_interactions(sim);
  
  // possibly do the RF contributions due to excluded atoms
  if(sim.nonbonded().RF_exclusion()){
    do_perturbed_RF_excluded_interactions(sim);
  }
  
}

/**
 * helper function to calculate perturbed forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
						 t_pairlist, t_innerloop>
::do_perturbed_interactions(t_simulation &sim,
			    typename t_pairlist::iterator it,
			    typename t_pairlist::iterator to)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  DEBUG(8, "Alpha LJ: " << m_nonbonded_interaction.alpha_lj() 
	<< " Alpha CRF: " << m_nonbonded_interaction.alpha_crf());

  for( ; it != to; ++it){    
    DEBUG(8, "perturbed pair: " << it.i() << " - " << *it);
    
    m_nonbonded_interaction.perturbed_interaction_inner_loop(sim, it.i(), *it);

  }

  // and long-range energy lambda-derivatives
  for(size_t i = 0; 
      i < m_nonbonded_interaction.pairlist().filter()
	.lambda_energies().lj_energy.size(); ++i){
    for(size_t j = 0; j < m_nonbonded_interaction.pairlist()
	  .filter().lambda_energies().lj_energy.size(); ++j){

      sim.system().lambda_energies().lj_energy[i][j] += 
	m_nonbonded_interaction.pairlist().filter().lambda_energies()
	.lj_energy[i][j];
      sim.system().lambda_energies().crf_energy[i][j] += 
	m_nonbonded_interaction.pairlist().filter().lambda_energies()
	.crf_energy[i][j];
    }
  }

  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
					   t_pairlist, t_innerloop>
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

      m_nonbonded_interaction.perturbed_one_four_interaction_inner_loop
	(sim, mit->second.sequence_number(), *it);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::do_perturbed_RF_excluded_interactions(t_simulation &sim)
{
  /*  
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
    m_nonbonded_interaction.rf_interaction(r,sim.topology().charge()(i) * sim.topology().charge()(i),
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
  */
}  

