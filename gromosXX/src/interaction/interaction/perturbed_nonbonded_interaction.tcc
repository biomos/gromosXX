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
template<typename t_simulation, typename t_pairlist, typename t_innerloop, typename t_nonbonded_interaction>
inline interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_pairlist, t_innerloop, t_nonbonded_interaction>
::Perturbed_Nonbonded_Interaction(t_simulation &sim, 
				  t_nonbonded_interaction &nonbonded_interaction)
  : Interaction<t_simulation>("Perturbed NonBonded"),
    m_nonbonded_interaction(nonbonded_interaction)
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_pairlist, 
	 typename t_innerloop, typename t_nonbonded_interaction>
inline interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_pairlist, t_innerloop, t_nonbonded_interaction>
::~Perturbed_Nonbonded_Interaction()
{
  DEBUG(4, "Perturbed_Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_pairlist, 
	 typename t_innerloop, typename t_nonbonded_interaction>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_pairlist, t_innerloop, t_nonbonded_interaction>
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

  // do the perturbed pairs
  do_perturbed_pair_interactions(sim);

}

/**
 * helper function to calculate perturbed forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_pairlist,
	 typename t_innerloop, typename t_nonbonded_interaction>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_pairlist, t_innerloop, t_nonbonded_interaction>
::do_perturbed_interactions(t_simulation &sim,
			    typename t_pairlist::iterator it,
			    typename t_pairlist::iterator to)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  for( ; it != to; ++it){    
    DEBUG(8, "ni: perturbed pair: " << it.i() << " - " << *it);
    
    m_nonbonded_interaction.perturbed_interaction_inner_loop(sim, it.i(), *it);

  }

  // and long-range energy lambda-derivatives
  DEBUG(7, "add long-range lambda-derivatives");

  for(size_t i = 0; 
      i < m_nonbonded_interaction.pairlist().filter()
	.lambda_energies().lj_energy.size(); ++i){
    for(size_t j = 0; j < m_nonbonded_interaction.pairlist()
	  .filter().lambda_energies().lj_energy.size(); ++j){

      assert(sim.system().lambda_energies().lj_energy.size() > i);
      assert(sim.system().lambda_energies().lj_energy[i].size() > j);
      assert(sim.system().lambda_energies().lj_energy.size() > j);
      assert(sim.system().lambda_energies().lj_energy[j].size() > i);
      
      sim.system().lambda_energies().lj_energy[i][j] += 
	m_nonbonded_interaction.pairlist().filter().lambda_energies()
	.lj_energy[i][j];
      sim.system().lambda_energies().crf_energy[i][j] += 
	m_nonbonded_interaction.pairlist().filter().lambda_energies()
	.crf_energy[i][j];
    }
  }

  DEBUG(7, "end of function perturbed nonbonded interaction");
  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_simulation, typename t_pairlist, 
	 typename t_innerloop, typename t_nonbonded_interaction>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_pairlist, t_innerloop, t_nonbonded_interaction>
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
template<typename t_simulation, typename t_pairlist, 
	 typename t_innerloop, typename t_nonbonded_interaction>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_pairlist, t_innerloop, t_nonbonded_interaction>
::do_perturbed_RF_excluded_interactions(t_simulation &sim)
{

  DEBUG(7, "\tcalculate perturbed excluded RF interactions");

  std::map<size_t, simulation::Perturbed_Atom>::const_iterator
    mit=sim.topology().perturbed_solute().atoms().begin(),
    mto=sim.topology().perturbed_solute().atoms().end();
  DEBUG(7, "\tSize of perturbed atoms " << sim.topology().perturbed_solute().atoms().size());
  
  for(; mit!=mto; ++mit){
    
    m_nonbonded_interaction.
      perturbed_RF_excluded_interaction_inner_loop(sim, mit);
    
  }
}

/**
 * calculate the interactions for the
 * PERTURBED PAIRS
 * (different interaction types in A and in B)
 */
template<typename t_simulation, typename t_pairlist,
	 typename t_innerloop, typename t_nonbonded_interaction>
inline void interaction::Perturbed_Nonbonded_Interaction<
  t_simulation, t_pairlist, t_innerloop, t_nonbonded_interaction>
::do_perturbed_pair_interactions(t_simulation &sim)
{

  std::vector<simulation::Perturbed_Atompair>::const_iterator
    it = sim.topology().perturbed_solute().atompairs().begin(),
    to = sim.topology().perturbed_solute().atompairs().end();
  
  math::Vec r, f, A_f, B_f;
  double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
    B_e_lj, B_e_crf, B_de_lj, B_de_crf;
  double e_lj, e_crf, de_lj, de_crf;
  lj_parameter_struct const *A_lj;
  lj_parameter_struct const *B_lj;
  double A_q, B_q;
  double alpha_lj, alpha_crf;
  
  const double l = sim.topology().lambda();
  DEBUG(7, "lambda: " << l);
  bool is_perturbed;
  
  for(; it != to; ++it){
    m_nonbonded_interaction.
      perturbed_pair_interaction_inner_loop(sim, it);
  }
  
}

