/**
 * @file simple_pairlist.tcc
 * simple pairlist implementation.
 */

template<typename t_simulation>
inline 
interaction::simple_pairlist<t_simulation>::iterator::iterator(
  t_pl_matrix &pl,
  int ii,
  int jj
)
  : 
  m_pairlist(pl),
  m_i(ii),
  m_j(jj)
{}

template<typename t_simulation>
inline 
void 
interaction::simple_pairlist<t_simulation>::iterator::operator++()
{
  ++j();
  if (j() == m_pairlist[i()].size()) {
    j() = 0;
    do {
      ++i();
      if (i() == m_pairlist.size()) 
        break;
    } while (!m_pairlist[i()].size());
  }
}

template<typename t_simulation>
inline 
bool
interaction::simple_pairlist<t_simulation>::iterator::operator==(
  typename simple_pairlist<t_simulation>::iterator it
)
{
  return (
    &pairlist() == &(it.pairlist()) &&
    i() == it.i() &&
    j() == it.j()
  );
}

template<typename t_simulation>
inline 
bool
interaction::simple_pairlist<t_simulation>::iterator::operator!=(
  typename simple_pairlist<t_simulation>::iterator it
)
{
  return !operator==(it);
}

template<typename t_simulation>
void interaction::simple_pairlist<t_simulation>
::update(t_simulation &sim)
{

  clear();
  size_t num_atoms = sim.topology().num_solute_atoms();
  resize(num_atoms);
  
  // now the a bit silly algorithm
  for(size_t i=0; i<num_atoms; ++i)
    for(size_t j=i+1; j<num_atoms; ++j){
      // check if not excluded
      if(i < sim.topology().solute().num_atoms())
	if (sim.topology().all_exclusion(i).count(j))
	  continue;

      (*this)[i].push_back(j);
    }
  
}

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

template<typename t_simulation>
void interaction::twin_range_pairlist_cg<t_simulation>
::update(t_simulation &sim)
{
  static math::VArray cg_cog;
  cg_cog.resize(sim.topology().num_chargegroups());

  short_range().clear();
  long_range().clear();

  size_t num_atoms = sim.topology().num_solute_atoms();

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

      sim.system().periodicity().nearest_image(cg_cog(cg1_index), cg_cog(cg2_index), p);
      d = dot(p, p);

      if (d > rcutl2)        // OUTSIDE
        continue;

      else if (d > rcutp2){  // LONGRANGE
	do_cg_interaction(cg1, cg2, long_range());
      }

      else {                 // SHORTRANGE
	if (unsigned(**cg2) < sim.topology().solute().num_atoms()){
	  // no exclusions... (at least cg2 is solvent)
	  do_cg_interaction(cg1, cg2, short_range());
	}
	else{
	  // exclusions!
	  do_cg_interaction_excl(sim, cg1, cg2, short_range());
	}
	
      } // ranges
    } // inter cg (cg2)
  } // cg1
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
      
      pl[*a1].push_back(*a2);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1
}


template<typename t_simulation>
inline 
typename interaction::simple_pairlist<t_simulation>::iterator 
interaction::simple_pairlist<t_simulation>::begin()
{
  for (unsigned int ii = 0; ii < size(); ii++)
    if ((*this)[ii].size())
      return iterator(*this, ii);
  return end();
}

template<typename t_simulation>
inline 
typename interaction::simple_pairlist<t_simulation>::iterator 
interaction::simple_pairlist<t_simulation>::end()
{
  return iterator(*this, size());
}

template<typename t_simulation>
std::ostream& 
interaction::operator<<(
  std::ostream &os, 
  class simple_pairlist<t_simulation>& pl
)
{
  os << "printing pairlist\n";
  os << "-----------------" << endl;

  typename simple_pairlist<t_simulation>::iterator it = pl.begin();
  while (it != pl.end()) {
    if (!it.j())
      os << endl << std::setw(6) << it.i() << " | " << flush;
    os << std::setw(6) << *it; 
    ++it;
  }
  os << std::endl;

  return os;
}
