/**
 * @file twinrange_pairlist_algorithm.tcc
 * implement methods of Twinrange_Pairlist_Algorithm.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include "../../debug.h"

template<typename t_simulation, typename t_filter>
interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, t_filter>
::Chargegroup_Range_Pairlist_Algorithm(basic_pairlist_type &pairlist,
				       basic_pairlist_type &perturbed_pl,
				       Nonbonded_Base &base)
  : Basic_Pairlist_Algorithm<t_simulation, t_filter>(pairlist, 
						     perturbed_pl, base) 
{
}

template<typename t_simulation, typename t_filter>
void interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, t_filter>
::update(t_simulation &sim)
{ 
  DEBUG(7, "pairlist update");
  
 
  // empty the pairlist
  m_pairlist.clear();
  m_pairlist.resize(sim.topology().num_atoms());

  // and the perturbed pairlist
  m_perturbed_pairlist.clear();
  m_perturbed_pairlist.resize(sim.topology().num_atoms());

  DEBUG(7, "pairlist resized");
  
  // prepare the filter
  m_filter.prepare(sim);
  
  DEBUG(7, "filter prepared");
  
  // loop over the chargegroups
  simulation::chargegroup_iterator cg1 = sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();
  
  for(int cg1_index=0; cg1 != cg_to; ++cg1, ++cg1_index) {
    // add intra cg (if not solvent...)
    DEBUG(7, "intra: " << cg1_index);
    if (unsigned(**cg1) < sim.topology().solute().num_atoms()){
      do_cg_interaction_intra(sim, cg1);
    }
  
    // inter chargegroup
    simulation::chargegroup_iterator cg2(*cg1+1);
    for(int cg2_index = cg1_index + 1; cg2 != cg_to; ++cg2, ++cg2_index) {

      DEBUG(7, "cg1: " << cg1_index << " - cg2: " << cg2_index);

      if (m_filter.range_chargegroup_pair(sim, cg1_index, cg2_index, cg1, cg2))
	continue;
      DEBUG(7, "it's in SHORT-range!");
      
      // SHORTRANGE
      if (unsigned(**cg2) < sim.topology().solute().num_atoms()){
	// exclusions!
	do_cg_interaction_excl(sim, cg1, cg2);
      }
      else{
	// no exclusions... (at least cg2 is solvent)
	do_cg_interaction(sim, cg1, cg2);
      }
	
    } // inter cg (cg2)
  } // cg1
  
  DEBUG(7, "pairlist done");
  DEBUG(7, "pairlist size: " << m_pairlist.size());
  
}

/**
 * inter cg, no exclusion
 */
template<typename t_simulation, typename t_filter>
inline void
interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, t_filter>
::do_cg_interaction(t_simulation const & sim,
		    simulation::chargegroup_iterator cg1,
		    simulation::chargegroup_iterator cg2)
{
  simulation::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      // check if one is perturbed
      if (m_filter.perturbed_atom(sim, *a1))
	{
	  // i perturbed (and maybe j)
	  m_perturbed_pairlist[*a1].push_back(*a2);
	  continue;
	}
      else if (m_filter.perturbed_atom(sim, *a2))
	{
	  m_perturbed_pairlist[*a2].push_back(*a1);
	  continue;
	}

#ifndef NDEBUG
      if (*a1 >= m_pairlist.size()){
	std::cout << "a1=" << *a1 << " a2=" << *a2
		  << " cg1=" << **cg1 << " cg2=" << **cg2
		  << " pairlist.size=" << m_pairlist.size() << std::endl;
      }
#endif
      assert(*a1 < m_pairlist.size());

      m_pairlist[*a1].push_back(*a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}


template<typename t_simulation, typename t_filter>
inline void
interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, t_filter>
::do_cg_interaction_excl(t_simulation &sim,
			 simulation::chargegroup_iterator cg1,
			 simulation::chargegroup_iterator cg2)
{
  simulation::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){
      DEBUG(7, "before exlusion check: " << *a1 << " -- " << *a2);
      
      // check it is not excluded
      if (m_filter.exclusion_solute_pair(sim, *a1, *a2)){
	DEBUG(7, "excluded" );
	continue;
      }

      // check if one is perturbed
      if (m_filter.perturbed_atom(sim, *a1))
	{
	  // i perturbed (and maybe j)
	  m_perturbed_pairlist[*a1].push_back(*a2);
	  continue;
	}
      else if (m_filter.perturbed_atom(sim, *a2))
	{
	  m_perturbed_pairlist[*a2].push_back(*a1);
	  continue;
	}
      
#ifndef NDEBUG
      if (*a1 >= m_pairlist.size()){
	std::cout << "a1=" << *a1 << " a2=" << *a2
		  << " cg1=" << **cg1 << " cg2=" << **cg2
		  << " pairlist.size=" << m_pairlist.size() << std::endl;
      }
#endif
      assert(*a1 < m_pairlist.size());

      m_pairlist[*a1].push_back(*a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_simulation, typename t_filter>
inline void
interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, t_filter>
::do_cg_interaction_intra(t_simulation &sim,
			  simulation::chargegroup_iterator cg1)
{
  simulation::Atom_Iterator a1 = cg1.begin(),
    a1_to = cg1.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2(*a1+1);
	a2 != a1_to; ++a2){

      // check it is not excluded
      if (m_filter.exclusion_solute_pair(sim, *a1, *a2))
	continue;
      // check if one is perturbed
      if (m_filter.perturbed_atom(sim, *a1))
	{
	  // i perturbed (and maybe j)
	  m_perturbed_pairlist[*a1].push_back(*a2);
	  continue;
	}
      else if (m_filter.perturbed_atom(sim, *a2))
	{
	  m_perturbed_pairlist[*a2].push_back(*a1);
	  continue;
	}

#ifndef NDEBUG
      if (*a1 >= m_pairlist.size()){
	std::cout << "a1=" << *a1 << " a2=" << *a2
		  << " cg1=" << **cg1
		  << " pairlist.size=" << m_pairlist.size() << std::endl;
      }
#endif
      assert(*a1 < m_pairlist.size());

      m_pairlist[*a1].push_back(*a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}
