/**
 * exclusion_filter.tcc
 * filter for exclusions.
 */

template<typename t_interaction_spec>
inline
interaction::Exclusion_Filter<t_interaction_spec>
::Exclusion_Filter()
  : interaction::Filter()
{
}

template<typename t_interaction_spec>
inline bool
interaction::Exclusion_Filter<t_interaction_spec>
::excluded_solute_pair(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       size_t const i, size_t const j)
{
  if (t_interaction_spec::do_exclusion){
    assert(i<j);

    std::set<int>::const_reverse_iterator
      e = topo.all_exclusion(i).rbegin(),
      e_to = topo.all_exclusion(i).rend();

    for( ; e != e_to; ++e){
      if (j > unsigned(*e)) break;
      if (j == unsigned(*e)){
	DEBUG(11, "\texcluded");
	return true;
      }
      
    }
    
    /*
    if (topo.all_exclusion(i).count(j)){
      DEBUG(11, "\texcluded");
      return true;
    }
    */
  }
  DEBUG(12, "\tnot excluded");
  return false;
}

template<typename t_interaction_spec>
inline bool
interaction::Exclusion_Filter<t_interaction_spec>
::inverse_excluded_solute_pair(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       size_t const i, size_t const j)
{
  if (t_interaction_spec::do_exclusion){
    assert(j<i);

    std::set<int>::const_reverse_iterator
      e = topo.all_exclusion(j).rbegin(),
      e_to = topo.all_exclusion(j).rend();

    for( ; e != e_to; ++e){
      if (i > *e) break;
      if (i == *e){
	DEBUG(11, "\texcluded");
	return true;
      }
    }

    /*
    if (topo.all_exclusion(j).count(i)){
      DEBUG(11, "\texcluded");
      return true;
    }
    */

  }
  DEBUG(12, "\tnot excluded");
  return false;
}
