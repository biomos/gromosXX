/**
 * @file forcefield.tcc
 * contains the inline functions for
 * forcefield.
 */

template<typename t_simulation>
inline interaction::forcefield<t_simulation>::forcefield()
{
}

template<typename t_simulation>
inline interaction::forcefield<t_simulation>::~forcefield()
{
  for(typename std::vector<interaction<t_simulation> *>::iterator
	it = m_interaction.begin(),
	to = m_interaction.end();
      it != to;
      ++it){
    delete *it;
  }
}

template<typename t_simulation>
inline void interaction::forcefield<t_simulation>
::add_interaction(interaction<t_simulation> *inter)
{
  m_interaction.push_back(inter);
}


template<typename t_simulation>
inline void interaction::forcefield<t_simulation>
::calculate_interactions(t_simulation &simu)
{
  for(typename std::vector<interaction<t_simulation> *>::iterator 
	it = m_interaction.begin(),
	to = m_interaction.end();
      it != to;
      ++it){
    (*it)->calculate_interactions(simu);
  }
}
