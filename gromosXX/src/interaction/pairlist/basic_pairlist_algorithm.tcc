/**
 * @File basic_pairlist_algorithm.tcc
 * template methods
 */

template<typename t_simulation, typename t_filter>
interaction::Basic_Pairlist_Algorithm<t_simulation, t_filter>
::Basic_Pairlist_Algorithm(basic_pairlist_type &pairlist,
			   basic_pairlist_type &perturbed_pl,
			   Nonbonded_Base &base)
  : m_pairlist(pairlist),
    m_perturbed_pairlist(perturbed_pl),
    m_filter(base)
{
}

template<typename t_simulation, typename t_filter>
t_filter & interaction::Basic_Pairlist_Algorithm<t_simulation, t_filter>
::filter()
{
  return m_filter;
}

template<typename t_simulation, typename t_filter>
void interaction::Basic_Pairlist_Algorithm<t_simulation, t_filter>
::update(t_simulation &sim)
{  
  const size_t num_atoms = sim.topology().num_atoms();
  const size_t num_solute_atoms = sim.topology().num_solute_atoms();
  size_t i, j;
 
  // empty the pairlist
  m_pairlist.clear();
  m_pairlist.resize(num_atoms);

  // and the perturbed pairlist
  m_perturbed_pairlist.clear();
  m_perturbed_pairlist.resize(num_atoms);
  
  // prepare the filter
  m_filter.prepare();
 
  // solute
  for(i=0; i<num_solute_atoms; ++i){
    for(j=i+1; j<num_solute_atoms; ++j){
      
      // check solute exclusion
      if (m_filter.exclusion_solute_pair(sim, i, j)) continue;
      // check if one is perturbed
      if (m_filter.perturbed_pair(sim, i))
	{
	  // i perturbed (and maybe j)
	  m_perturbed_pairlist[i].push_back(j);
	  continue;
	}
      else if (m_filter.perturbed_atom(sim, j))
	{
	  m_perturbed_parlist[j].push_back(i);
	  continue;
	}
      
      m_pairlist[i].push_back(j);
    }
    // solute - solvent
    for( ; j < num_atoms; ++j){

      // check if one is perturbed
      if (m_filter.perturbed_pair(sim, i))
	{
	  // i perturbed (and maybe j)
	  m_perturbed_pairlist[i].push_back(j);
	  continue;
	}
      else if (m_filter.perturbed_atom(sim, j))
	{
	  m_perturbed_pairlist[j].push_back(i);
	  continue;
	}

      m_pairlist[i].push_back(j);
    }
  }
  
  // solvent
  for( ; i<num_atoms; ++i){
    for(j=i+1; j<num_atoms; ++j){
      
      // check solvent exclusions
      m_filter.exclusion_solvent_pair(sim, i, j);
      
      m_pairlist[i].push_back(j);
    }
  }
}

  
