/**
 * @file forcefield.tcc
 * contains the inline functions for
 * forcefield.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE forcefield

#include "../../debug.h"

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
::calculate_interactions(t_simulation &sim)
{
  DEBUG(5, "forcefield: calculate interaction");

  sim.system().force() = 0.0;
  sim.system().energies().zero();

  for(typename Forcefield<t_simulation>::iterator it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(7, "interaction: " << (*it)->name);
    (*it)->calculate_interactions(sim);
  }
}
