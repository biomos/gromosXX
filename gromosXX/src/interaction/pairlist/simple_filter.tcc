/**
 * @file simple_filter.tcc
 * simple filter methods.
 */

template<typename t_simulation>
inline bool
interaction::Simple_Filter<t_simulation>::solute_pair(t_simulation const &sim,
						      size_t const i,
						      size_t const j)
{
  return filter_excluded(sim, i, j);
}

template<typename t_simulation>
inline bool
interaction::Simple_Filter<t_simulation>::solvent_pair(t_simulation const &sim,
						       size_t const i,
						       size_t const j)
{
  return filter_solvent_excluded(sim, i, j);
}

template<typename t_simulation>
inline bool
interaction::Simple_Filter<t_simulation>
::filter_excluded(t_simulation const & sim, size_t const i, size_t const j)
{
  // check it is not excluded
  if (sim.topology().all_exclusion(i).count(j))
    return true;
  return false;
}

template<typename t_simulation>
inline bool
interaction::Simple_Filter<t_simulation>
::filter_solvent_excluded(t_simulation const & sim,
			  size_t i, size_t j)
{
  size_t s = 0;
  for( ; i - sim.topology().num_solvent_atoms[s] > 0;
       i -= sim.topology().num_solvent_atoms[s++]){
  }
  i /= sim.topology().solvent(s-1).num_atoms();
  
  s = 0;
  for( ; j - sim.topology().num_solvent_atoms[s] > 0;
       j -= sim.topology().num_solvent_atoms[s++]){
  }
  j /= sim.topology().solvent(s-1).num_atoms();
  
  if (i == j) return true;

  return false;
    
}
