/**
 * @file twin_range_pairlist.tcc
 * twin range
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include "../../debug.h"

template<typename t_simulation>
void interaction::twin_range_pairlist<t_simulation>
::update(t_simulation &sim)
{
  short_range().clear();
  long_range().clear();

  size_t num_atoms = sim.topology().num_atoms();
  DEBUG(10, "pairlist size:" << num_atoms << " solute: " << sim.topology().num_solute_atoms() << " solvent: " << sim.topology().num_solvent_atoms());
  

  short_range().resize(num_atoms);
  long_range().resize(num_atoms);

  math::VArray &pos = sim.system().pos();

  double d = 0;

  // square of the cutoff...
  double rcutp2 = sim.nonbonded().cutoff_short();
  rcutp2 *= rcutp2;
  
  double rcutl2 = sim.nonbonded().cutoff_long();
  rcutl2 *= rcutl2;

  math::Vec p;
  
  for(size_t i=0; i<num_atoms; ++i) {
    for(size_t j=i+1; j<num_atoms; ++j) {

      // check it is not excluded
      if (sim.topology().all_exclusion(i).count(j))
        continue;

      sim.system().periodicity().nearest_image(pos(i), pos(j), p);
      d = dot(p, p);

      if (d > rcutl2) 
        continue;
      else if (d > rcutp2)
        long_range()[i].push_back(j);
      else {
        short_range()[i].push_back(j);
      }
    }
  }
}
