/**
 * @file basic_pairlist_algorithm.tcc
 * template methods
 */

template<typename t_simulation, typename t_filter>
interaction::Basic_Pairlist_Algorithm<t_simulation, t_filter>
::Basic_Pairlist_Algorithm(std::vector<std::vector<unsigned int> > &pairlist)
  : m_pairlist(pairlist)
{
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
 
  // solute
  for(i=0; i<num_solute_atoms; ++i){
    for(j=i+1; j<num_solute_atoms; ++j){
      
      // check solute exclusion
      
      m_pairlist[i].push_back(j);
    }
    // solute - solvent
    for( ; j < num_atoms; ++j){
      m_pairlist[i].push_back(j);
    }
  }
  
  // solvent
  for( ; i<num_atoms; ++i){
    for(j=i+1; j<num_atoms; ++j){
      
      // check solvent exclusions
      
      m_pairlist[i].push_back(j);
    }
  }
}

  
