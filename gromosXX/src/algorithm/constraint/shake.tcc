/**
 * @file shake.tcc
 * contains the template methods for
 * the class shake.
 */

/**
 * Constructor.
 */
template<typename t_simulation>
algorithm::shake<t_simulation>
::shake()
  : tolerance(0.0001),
    max_iterations(5)
{
}

/**
 * shake solvent.
 */
template<typename t_simulation>
int algorithm::shake<t_simulation>
::solvent(typename simulation_type::topology_type &topo,
	  typename simulation_type::system_type &sys,
	  double dt)
{
  // the first atom of a solvent
  size_t first = topo.num_solute_atoms();

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int tot_iterations = 0;

  // for all solvents
  for(size_t i=0; i<topo.num_solvents(); ++i){

    // loop over the molecules
    for(size_t nm=0; nm<topo.num_solvent_molecules(i);
	++nm, first+=topo.solvents(i).num_atoms()){

      skip_now.assign(topo.solvents(i).num_atoms(), false);
      skip_next.assign(topo.solvents(i).num_atoms(), true);

      int num_iterations = 0;
      bool convergence = false;
      while(!convergence){
	convergence = true;
	
	// and constraints
	for(std::vector<simulation::solvent::solventconstraint_struct>
	      ::const_iterator
	      it = topo.solvents(i).constraints().begin(),
	      to = topo.solvents(i).constraints().end();
	    it != to;
	    ++it){
	
	  // check whether we can skip this constraint
	  if (skip_now[it->i] && skip_now[it->j]) continue;

	  // std::cout << "mol: " << nm << " i: " << it->i << " j: " << it->j
	  // << " first: " << first << std::endl;

	  // the position
	  math::Vec &pos_i = sys.pos()(first+it->i);
	  math::Vec &pos_j = sys.pos()(first+it->j);
	
	  math::Vec r = pos_i - pos_j;
	  double dist2 = dot(r, r);
	
	  double constr_length2 = it->b0 * it->b0;
	  double diff = constr_length2 - dist2;

	  // std::cout << "constr: " << constr_length2 << " dist2: " << dist2 << std::endl;
	  
	  if(fabs(diff) >= constr_length2 * tolerance * 2.0){
	    // we have to shake

	    // the reference position
	    const math::Vec &ref_i = sys.old_pos()(first+it->i);
	    const math::Vec &ref_j = sys.old_pos()(first+it->j);

	    math::Vec ref_r = ref_i - ref_j;
	    double sp = dot(ref_r, r);
	  
	    if(sp < constr_length2 * math::epsilon){
	      io::messages.add("SHAKE error. vectors orthogonal",
			       "shake::solvent",
			       io::message::critical);
	      throw std::runtime_error("SHAKE failure in solvent");
	    }
	  
	    // lagrange multiplier
	    double lambda = diff / (sp * 2 *
				    (1.0 / topo.mass()(first+it->i) +
				     1.0 / topo.mass()(first+it->j) ));
	    // update positions
	    ref_r *= lambda;
	    pos_i += ref_r / topo.mass()(first+it->i);
	    pos_j -= ref_r / topo.mass()(first+it->j);
	  
	    convergence = false;

	    // consider atoms in the next step
	    skip_next[it->i] = false;
	    skip_next[it->j] = false;
	  
	  } // we have to shake
	} // constraints
      
	skip_now = skip_next;
	skip_next.assign(topo.solvents(i).num_atoms(), true);
      
	// std::cout << num_iterations+1 << std::endl;
	if(++num_iterations > max_iterations){
	  io::messages.add("SHAKE error. too many iterations",
			  "shake::solvent",
			  io::message::critical);
	  throw std::runtime_error("SHAKE failure in solvent");
	}
	
      } // while(!convergence)
      
      tot_iterations += num_iterations;
      
    } // molecules
    
  } // solvents

  // VELOCITY CORRECTION
  sys.vel() = (sys.pos() - sys.old_pos()) / dt;

  return tot_iterations;
  
} // shake solvent

