/**
 * @file twin_range_pairlist_cg.tcc
 * twin range with chargegroups
 */


template<typename t_simulation>
void interaction::twin_range_pairlist_cg<t_simulation>
::update(t_simulation &sim)
{
  short_range().clear();
  long_range().clear();

  static math::VArray cg_cog;
  cg_cog.resize(sim.topology().num_chargegroups());

  size_t num_atoms = sim.topology().num_atoms();

  DEBUG(10, "pairlist size:" << num_atoms << " solute: " << sim.topology().num_solute_atoms() << " solvent: " << sim.topology().num_solvent_atoms());
  
  short_range().resize(num_atoms);
  long_range().resize(num_atoms);

  math::VArray &pos = sim.system().pos();

  // calculate all center of geometries
  simulation::chargegroup_iterator cg1 = sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();

  size_t i;
  // math::Vec v;
  
  for(i=0; i < sim.topology().num_solute_chargegroups(); ++cg1, ++i){
    cg1.cog(pos, cg_cog(i));
    // cg_cog(i) = v;
  }
  for( ; cg1 != cg_to; ++cg1, ++i){
    cg_cog(i) = pos(**cg1);
  }

  double d = 0;

  // square of the cutoff...
  double rcutp2 = sim.nonbonded().cutoff_short();
  rcutp2 *= rcutp2;
  
  double rcutl2 = sim.nonbonded().cutoff_long();
  rcutl2 *= rcutl2;

  math::Vec cog1, cog2;
  math::Vec p;
  
  // loop over the chargegroups
  cg1 = sim.topology().chargegroup_begin();
  
  for(int cg1_index=0; cg1 != cg_to; ++cg1, ++cg1_index) {
    // add intra cg (if not solvent...)
    if (unsigned(**cg1) < sim.topology().solute().num_atoms()){
      do_cg_interaction_intra(sim, cg1, short_range());
    }
    
    // inter chargegroup
    simulation::chargegroup_iterator cg2(*cg1+1);
    for(int cg2_index = cg1_index + 1; cg2 != cg_to; ++cg2, ++cg2_index) {

      sim.system().periodicity().
	nearest_image(cg_cog(cg1_index), cg_cog(cg2_index), p);
      d = dot(p, p);

      if (d > rcutl2)        // OUTSIDE
        continue;

      else if (d > rcutp2){  // LONGRANGE
	do_cg_interaction(cg1, cg2, long_range());
      }

      else {                 // SHORTRANGE
	if (unsigned(**cg2) < sim.topology().solute().num_atoms()){
	  // exclusions!
	  do_cg_interaction_excl(sim, cg1, cg2, short_range());
	}
	else{
	  // no exclusions... (at least cg2 is solvent)
	  do_cg_interaction(cg1, cg2, short_range());
	}
	
      } // ranges
    } // inter cg (cg2)
  } // cg1

  DEBUG(7, "pairlist done");
  
}

template<typename t_simulation>
static void do_cg_interaction(simulation::chargegroup_iterator cg1,
			      simulation::chargegroup_iterator cg2,
			      interaction::simple_pairlist<t_simulation> &pl)
{
  simulation::chargegroup_iterator::atom_iterator a1 = cg1.begin(),
    a1_to = cg1.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::chargegroup_iterator::atom_iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      if (*a1 >= pl.size()){
	std::cout << "a1=" << *a1 << " a2=" << *a2
		  << " cg1=" << **cg1 << " cg2=" << **cg2
		  << " pl.size=" << pl.size() << std::endl;
      }

      assert(*a1 < pl.size());
      pl[*a1].push_back(*a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_simulation>
static void do_cg_interaction_excl(t_simulation &sim,
				   simulation::chargegroup_iterator cg1,
				   simulation::chargegroup_iterator cg2,
				   interaction::simple_pairlist<t_simulation> &pl)
{
  simulation::chargegroup_iterator::atom_iterator a1 = cg1.begin(),
    a1_to = cg1.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::chargegroup_iterator::atom_iterator
	  a2 = cg2.begin(),
	  a2_to = cg2.end();
	a2 != a2_to; ++a2){

      // check it is not excluded
      if(*a2 < sim.topology().solute().num_atoms())
	if (sim.topology().all_exclusion(*a1).count(*a2))
	  continue;

      assert(*a1 < pl.size());
      pl[*a1].push_back(*a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}

template<typename t_simulation>
static void do_cg_interaction_intra(t_simulation &sim,
				   simulation::chargegroup_iterator cg1,
				   interaction::simple_pairlist<t_simulation> &pl)
{
  simulation::chargegroup_iterator::atom_iterator a1 = cg1.begin(),
    a1_to = cg1.end();
  
  for( ; a1 != a1_to; ++a1){
    for(simulation::chargegroup_iterator::atom_iterator
	  a2(*a1+1);
	a2 != a1_to; ++a2){

      // check it is not excluded
      if (sim.topology().all_exclusion(*a1).count(*a2))
	continue;
      
      assert(*a1 < pl.size());
      pl[*a1].push_back(*a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}
