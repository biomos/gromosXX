/**
 * @file forcefield.tcc
 * contains the inline functions for
 * forcefield.
 */

template<typename t_simulation>
inline interaction::Forcefield<t_simulation>::Forcefield()
  : std::vector<Interaction<t_simulation> *>()
{
}

template<typename t_simulation>
inline interaction::Forcefield<t_simulation>::~Forcefield()
{
  for(typename Forcefield<t_simulation>::iterator it = begin(), to = end();
      it != to;
      ++it){
    delete *it;
  }
}

template<typename t_simulation>
inline void interaction::Forcefield<t_simulation>
::calculate_interactions(t_simulation &simu)
{

  simu.system().force() = 0.0;

  for(typename Forcefield<t_simulation>::iterator it = begin(), to = end();
      it != to;
      ++it){
    (*it)->calculate_interactions(simu);
  }
}
