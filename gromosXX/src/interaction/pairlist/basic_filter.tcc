/**
 * @file basic_filter.tcc
 * basic filter methods.
 */

template<typename t_simulation, typename t_base>
inline
interaction::Basic_Filter<t_simulation, t_base>
::Basic_Filter(t_base &base)
  : m_base(base)
{
}

template<typename t_simulation, typename t_base>
inline void
interaction::Basic_Filter<t_simulation, t_base>
::prepare(t_simulation &sim)
{
}

template<typename t_simulation, typename t_base>
inline bool
interaction::Basic_Filter<t_simulation, t_base>
::exclusion_solute_pair(t_simulation const &sim,
			size_t const i,
			size_t const j)
{
  // check it is not excluded
  if (sim.topology().all_exclusion(i).count(j))
    return true;
  return false;
}

template<typename t_simulation, typename t_base>
inline bool
interaction::Basic_Filter<t_simulation, t_base>
::exclusion_solvent_pair(t_simulation const &sim,
			 size_t const i,
			 size_t const j)
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
