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

template<typename t_simulation, typename t_interaction_spec>
inline interaction::Forcefield<t_simulation, t_interaction_spec>::Forcefield()
  : std::vector<Interaction<t_simulation, t_interaction_spec> *>()
{
}

template<typename t_simulation, typename t_interaction_spec>
inline interaction::Forcefield<t_simulation, t_interaction_spec>::~Forcefield()
{
  for(typename Forcefield<t_simulation, t_interaction_spec>::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    delete *it;
  }
}

template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Forcefield<t_simulation, t_interaction_spec>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(5, "forcefield: calculate interaction");

  sim.system().force() = 0.0;
  DEBUG(15, "zero energies");
  sim.system().energies().zero();
  DEBUG(15, "zero lambda energies");
  sim.system().lambda_energies().zero();
  sim.system().virial() = 0.0;

  // prepare for the virial
  if (t_interaction_spec::do_virial != interaction::no_virial){
    if(sim.pressure_calculation())
      sim.calculate_mol_com();
  }
  
  for(typename Forcefield<t_simulation, t_interaction_spec>::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(7, "interaction: " << (*it)->name);
    (*it)->calculate_interactions(sim);
  }

  // prefactor to the virial
  // done before the pressure calculation...
  /*
  if (t_interaction_spec::do_virial != interaction::no_virial){
    for(size_t i=0; i<3; ++i)
      for(size_t j=0; j<3; ++j)
	sim.system().virial()(i,j) =
	  -0.5 * sim.system().virial()(i,j);
  }
  */

}
