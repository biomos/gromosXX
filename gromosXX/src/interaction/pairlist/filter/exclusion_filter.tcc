/**
 * exclusion_filter.tcc
 * filter for exclusions.
 */

template<typename t_simulation, typename t_nonbonded_spec>
inline
interaction::Exclusion_Filter<t_simulation, t_nonbonded_spec>
::Exclusion_Filter()
  : interaction::Filter<t_simulation, t_nonbonded_spec>()
{
}

template<typename t_simulation, typename t_nonbonded_spec>
inline bool
interaction::Exclusion_Filter<t_simulation, t_nonbonded_spec>
::excluded_solute_pair(t_simulation const &sim,
		       size_t const i, size_t const j)
{
  if (t_nonbonded_spec::do_exclusion){
    assert(i<j);
    if (sim.topology().all_exclusion(i).count(j)){
      DEBUG(11, "\texcluded");
      return true;
    }
  }
  DEBUG(11, "\tnot excluded");
  return false;
}
