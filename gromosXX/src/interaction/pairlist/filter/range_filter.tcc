/**
 * @file range_filter.tcc
 * methods of Range_Filter
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include "../../../debug.h"


template<typename t_simulation, typename t_nonbonded_spec>
inline
interaction::Range_Filter<t_simulation, t_nonbonded_spec>
::Range_Filter()
  : interaction::Filter<t_simulation, t_nonbonded_spec>()
{
}

template<typename t_simulation, typename t_nonbonded_spec>
inline void
interaction::Range_Filter<t_simulation, t_nonbonded_spec>
::set_cutoff(double const cutoff_short, double const cutoff_long)
{
  m_cutoff_long = cutoff_long;
  m_cutoff_short = cutoff_short;
  m_cutoff_short_2 = cutoff_short * cutoff_short;
  m_cutoff_long_2  = cutoff_long * cutoff_long;
}

template<typename t_simulation, typename t_nonbonded_spec>
inline void
interaction::Range_Filter<t_simulation, t_nonbonded_spec>
::prepare_cog(t_simulation & sim)
{
  DEBUG(11, "Range_Filter::prepare_cog");
  
  // first put the chargegroups into the box
  sim.system().periodicity().put_chargegroups_into_box(sim);

  // calculate cg cog's
  m_cg_cog.resize(sim.topology().num_chargegroups());
  math::VArray const &pos = sim.system().pos();

  // calculate all center of geometries
  simulation::chargegroup_iterator
    cg1 =   sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();

  size_t i;

  // solute
  for(i=0; i < sim.topology().num_solute_chargegroups(); ++cg1, ++i){
    cg1.cog(pos, m_cg_cog(i));
  }
  // solvent
  for( ; cg1 != cg_to; ++cg1, ++i){
    m_cg_cog(i) = pos(**cg1);
  }  
}

template<typename t_simulation, typename t_nonbonded_spec>
inline void
interaction::Range_Filter<t_simulation, t_nonbonded_spec>
::grid_cog(t_simulation const & sim,
	   Chargegroup_Grid<t_simulation> & grid)
{
  DEBUG(8, "box:\n\t" << sim.system().periodicity().box()(0)(0)
	<< "\t" << sim.system().periodicity().box()(0)(1)
	<< "\t" << sim.system().periodicity().box()(0)(2)
	<< "\n\t" << sim.system().periodicity().box()(1)(0)
	<< "\t" << sim.system().periodicity().box()(1)(1)
	<< "\t" << sim.system().periodicity().box()(1)(2)
	<< "\n\t" << sim.system().periodicity().box()(2)(0)
	<< "\t" << sim.system().periodicity().box()(2)(1)
	<< "\t" << sim.system().periodicity().box()(2)(2));
	
  for(size_t i = 0, i_to = m_cg_cog.size(); i != i_to; ++i)
    grid.add(m_cg_cog(i), i);
  
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline bool
interaction::Range_Filter<t_simulation, t_nonbonded_spec>
::range_chargegroup_pair(t_simulation & sim,
			 t_nonbonded_interaction & nonbonded_interaction,
			 size_t const i, size_t const j,
			 simulation::chargegroup_iterator const & it_i,
			 simulation::chargegroup_iterator const & it_j)
{
  DEBUG(11, "Range_Filter::range_chargegroup_pair " << i << " - " << j);
  
  math::Vec p;
  
  assert(unsigned(m_cg_cog.size()) > i &&
	 unsigned(m_cg_cog.size()) > j);
  
  sim.system().periodicity().
    nearest_image(m_cg_cog(i), m_cg_cog(j), p);
 
  // the distance
  const double d = dot(p, p);

  DEBUG(11, "\tdistance: " << d);

  if (d > m_cutoff_long_2){        // OUTSIDE: filter
    DEBUG(11, "cg pair " << i << " - " << j << " outside range");
    return true;
  }
  
  if (d < m_cutoff_short_2){       // SHORTRANGE: no filter
    DEBUG(11, "cg pair " << i << " - " << j << " short range");
    return false;
  }

  DEBUG(11, "cg pair " << i << " - " << j << " long range");  
    
  // LONGRANGE: interactions and filter

  simulation::Atom_Iterator a1 = it_i.begin(),
    a1_to = it_i.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2 = it_j.begin(),
	  a2_to = it_j.end();
	a2 != a2_to; ++a2){

      // the interactions
      nonbonded_interaction.add_longrange_pair(sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1

  return true;
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline bool
interaction::Range_Filter<t_simulation, t_nonbonded_spec>
::range_chargegroup_pair(t_simulation & sim,
			 t_nonbonded_interaction & nonbonded_interaction,
			 size_t const i, size_t const j,
			 simulation::chargegroup_iterator const & it_i,
			 simulation::chargegroup_iterator const & it_j,
			 typename math::Boundary_Implementation
			 <t_simulation::system_type::boundary_type>
			 ::shift_struct const & shift)
{
  DEBUG(11, "Range_Filter::range_chargegroup_pair (shift) " << i << " - " << j);
  
  math::Vec p;
  
  assert(unsigned(m_cg_cog.size()) > i &&
	 unsigned(m_cg_cog.size()) > j);
  
  p = m_cg_cog(i) + shift.pos - m_cg_cog(j);
  
  // the distance
  const double d = dot(p, p);
  DEBUG(11, "\tdistance: " << d);

  if (d > m_cutoff_long_2){        // OUTSIDE: filter
    DEBUG(11, "cg pair " << i << " - " << j << " outside range");
    return true;
  }
  
  if (d < m_cutoff_short_2){       // SHORTRANGE: no filter
    DEBUG(11, "cg pair " << i << " - " << j << " short range");
    return false;
  }

  DEBUG(11, "cg pair " << i << " - " << j << " long range");  
  // LONGRANGE: interactions and filter

  simulation::Atom_Iterator a1 = it_i.begin(),
    a1_to = it_i.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::Atom_Iterator
	  a2 = it_j.begin(),
	  a2_to = it_j.end();
	a2 != a2_to; ++a2){

      // the interactions
      nonbonded_interaction.add_longrange_pair(sim, *a1, *a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1

  return true;
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline bool
interaction::Range_Filter<t_simulation, t_nonbonded_spec>
::range_atom_pair(t_simulation & sim,
		  t_nonbonded_interaction &nonbonded_interaction,
		  size_t const i, size_t const j)
{

  DEBUG(11, "Range_Filter::range_atom_pair " << i << " - " << j);
  
  math::Vec p;
  
  sim.system().periodicity().
    nearest_image(sim.system().pos()(i), sim.system().pos()(j), p);
 
  // the distance
  const double d = dot(p, p);

  if (d > m_cutoff_long_2){        // OUTSIDE: filter
    DEBUG(11, "atom pair " << i << " - " << j << " outside range");
    return true;
  }
  
  if (d < m_cutoff_short_2){       // SHORTRANGE: no filter
    DEBUG(11, "atom pair " << i << " - " << j << " short range");
    return false;
  }

  DEBUG(11, "atom pair " << i << " - " << j << " long range");  
    
  // LONGRANGE: interactions and filter

  // the interactions
  nonbonded_interaction.add_longrange_pair(sim, i, j);

  return true;
}
