/**
 * exclusion_filter.tcc
 * filter for exclusions.
 */

template<typename t_nonbonded_spec>
inline
interaction::Exclusion_Filter<t_nonbonded_spec>
::Exclusion_Filter()
  : interaction::Filter<t_nonbonded_spec>()
{
}

template<typename t_nonbonded_spec>
inline bool
interaction::Exclusion_Filter<t_nonbonded_spec>
::excluded_solute_pair(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       size_t const i, size_t const j)
{
  if (t_nonbonded_spec::do_exclusion){
    assert(i<j);
    if (topo.all_exclusion(i).count(j)){
      DEBUG(11, "\texcluded");
      return true;
    }
  }
  DEBUG(12, "\tnot excluded");
  return false;
}

template<typename t_nonbonded_spec>
inline bool
interaction::Exclusion_Filter<t_nonbonded_spec>
::inverse_excluded_solute_pair(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       size_t const i, size_t const j)
{
  if (t_nonbonded_spec::do_exclusion){
    assert(j<i);
    if (topo.all_exclusion(j).count(i)){
      DEBUG(11, "\texcluded");
      return true;
    }
  }
  DEBUG(12, "\tnot excluded");
  return false;
}
