/**
 * @file twinrange_chargegroup_filter.tcc
 * methods of Twinrange_Chargegroup_Filter
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include "../../debug.h"


template<typename t_simulation, typename t_base, typename t_innerloop, typename t_basic_filter>
inline
interaction::Twinrange_Chargegroup_Filter<t_simulation, t_base, t_innerloop, t_basic_filter>
::Twinrange_Chargegroup_Filter(t_base & base)
  : Twinrange_Filter<t_simulation, t_base, t_innerloop, t_basic_filter>(base)
{
}

template<typename t_simulation, typename t_base, typename t_innerloop, typename t_basic_filter>
inline void
interaction::Twinrange_Chargegroup_Filter<t_simulation, t_base, t_innerloop, t_basic_filter>
::prepare(t_simulation &sim)
{
  // call parent
  Twinrange_Filter<t_simulation, t_base, t_innerloop, 
    t_basic_filter>
    ::prepare(sim);
  
  // calculate cg cog's
  m_cg_cog.resize(sim.topology().num_chargegroups());
  math::VArray const &pos = sim.system().pos();

  // calculate all center of geometries
  simulation::chargegroup_iterator cg1 = sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();

  size_t i;
  
  for(i=0; i < sim.topology().num_solute_chargegroups(); ++cg1, ++i){
    cg1.cog(pos, m_cg_cog(i));
  }
  for( ; cg1 != cg_to; ++cg1, ++i){
    m_cg_cog(i) = pos(**cg1);
  }  
}

template<typename t_simulation, typename t_base, typename t_innerloop, typename t_basic_filter>
inline bool
interaction::Twinrange_Chargegroup_Filter<t_simulation, t_base, t_innerloop, t_basic_filter>
::range_chargegroup_pair(t_simulation &sim,
			 size_t const i, size_t const j,
			 simulation::chargegroup_iterator const & it_i,
			 simulation::chargegroup_iterator const & it_j)
{

  math::Vec p;
  
  sim.system().periodicity().
    nearest_image(m_cg_cog(i), m_cg_cog(j), p);
  const double d = dot(p, p);

  if (d > m_cutoff_long_2){        // OUTSIDE: filter
    DEBUG(7, "it's outside range");
    return true;
  }
  
  if (d < m_cutoff_short_2){       // SHORTRANGE: no filter
    DEBUG(7, "it's short range");
    return false;
  }
  
  DEBUG(7, "it's LONG-range");
  
  // LONGRANGE: interactions and filter

  simulation::Atom_Iterator a1 = it_i.begin(),
    a1_to = it_i.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2 = it_j.begin(),
	  a2_to = it_j.end();
	a2 != a2_to; ++a2){

      // the interactions


      if (perturbed_atom(sim, *a1))
	{
	  // i perturbed (and maybe j)
	  perturbed_interaction_inner_loop(sim, *a1, *a2);
	  continue;
	}
      else if (perturbed_atom(sim, *a2))
	{
	  perturbed_interaction_inner_loop(sim, *a2, *a1);
	  continue;
	}

      DEBUG(7, "inner loop: " << *a1 << " - " << *a2);      
      interaction_inner_loop(sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1

  return true;
}
