/**
 * @file twin_range_pairlist.tcc
 * twin range
 */

template<typename t_simulation>
void interaction::twin_range_pairlist<t_simulation>
::update(t_simulation &sim)
{
  short_range().clear();
  long_range().clear();

  size_t num_atoms = sim.topology().num_solute_atoms();

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

      sim.system().periodicity().nearest_image(pos(i), pos(j), p);
      d = dot(p, p);

      if (d > rcutl2) 
        continue;
      else if (d > rcutp2)
        long_range()[i].push_back(j);
      else {
	// check it is not excluded
	if(i < sim.topology().solute().num_atoms())
	  if (sim.topology().all_exclusion(i).count(j))
	    continue;

        short_range()[i].push_back(j);
      }
    }
  }
}
