/**
 * @file shake.tcc
 * contains the template methods for
 * the class Shake.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraint

#include "../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation>
algorithm::Shake<t_simulation>
::Shake(double const tolerance, int const max_iterations)
  : m_tolerance(tolerance),
    max_iterations(max_iterations)
{
}

template<typename t_simulation>
void algorithm::Shake<t_simulation>
::tolerance(double const tol)
{
  m_tolerance = tol;
}

/**
 * shake solute
 */
template<typename t_simulation>
int algorithm::Shake<t_simulation>
::solute(typename simulation_type::topology_type &topo,
	 typename simulation_type::system_type &sys,
	 double dt)
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...
  
  // check whether we shake
  if (!topo.solute().distance_constraints().size()) return 0;

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  int first = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);
  
  bool convergence = false;
  while(!convergence){
    convergence = _shake(topo, sys, first, skip_now, skip_next,
			 topo.solute().distance_constraints());

    // std::cout << num_iterations+1 << std::endl;
    if(++num_iterations > max_iterations){
      io::messages.add("SHAKE error. too many iterations",
		       "Shake::solute",
		       io::message::critical);
      throw std::runtime_error("SHAKE failure in solute");
    }

  } // convergence?

  // shaken velocity
  sys.vel() = (sys.pos() - sys.old_pos()) / dt;

  return num_iterations;

} // solute


/**
 * shake solvent.
 */
template<typename t_simulation>
int algorithm::Shake<t_simulation>
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
	++nm, first+=topo.solvent(i).num_atoms()){

      skip_now.assign(topo.solvent(i).num_atoms(), false);
      skip_next.assign(topo.solvent(i).num_atoms(), true);

      int num_iterations = 0;
      bool convergence = false;
      while(!convergence){
	convergence = _shake(topo, sys, first, skip_now, skip_next,
			     topo.solvent(i).distance_constraints());
	
	// std::cout << num_iterations+1 << std::endl;
	if(++num_iterations > max_iterations){
	  io::messages.add("SHAKE error. too many iterations",
			   "Shake::solvent",
			   io::message::critical);
	  throw std::runtime_error("SHAKE failure in solvent");
	}


      } // while(!convergence)
      
      tot_iterations += num_iterations;
      
    } // molecules
    
  } // solvents

  // shaken velocity
  sys.vel() = (sys.pos() - sys.old_pos()) / dt;

  return tot_iterations;
  
} // shake solvent

/**
 * do one iteration
 */      
template<typename t_simulation>
bool algorithm::Shake<t_simulation>
::_shake(typename simulation_type::topology_type const &topo,
	 typename simulation_type::system_type &sys,
	 int const first,
	 std::vector<bool> &skip_now,
	 std::vector<bool> &skip_next,
	 std::vector<simulation::compound::distance_constraint_struct>
	 & constr)
{
  bool convergence = true;

  // and constraints
  for(std::vector<simulation::compound::distance_constraint_struct>
	::const_iterator
	it = constr.begin(),
	to = constr.end();
      it != to;
      ++it){
	
    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j]) continue;

    DEBUG(10, "i: " << it->i << " j: " << it->j << " first: " << first);

    // the position
    math::Vec &pos_i = sys.pos()(first+it->i);
    math::Vec &pos_j = sys.pos()(first+it->j);

    DEBUG(10, "i: " << pos_i << "\nj: " << pos_j);
	
    math::Vec r = pos_i - pos_j;
    double dist2 = dot(r, r);
	
    double constr_length2 = it->b0 * it->b0;
    double diff = constr_length2 - dist2;

    DEBUG(15, "constr: " << constr_length2 << " dist2: " << dist2);
	  
    if(fabs(diff) >= constr_length2 * m_tolerance * 2.0){
      // we have to shake
      
      // the reference position
      const math::Vec &ref_i = sys.old_pos()(first+it->i);
      const math::Vec &ref_j = sys.old_pos()(first+it->j);

      math::Vec ref_r = ref_i - ref_j;
      double sp = dot(ref_r, r);
	  
      if(sp < constr_length2 * math::epsilon){
	io::messages.add("SHAKE error. vectors orthogonal",
			 "Shake::???",
			 io::message::critical);
	DEBUG(5, "ref i " << ref_i << " ref j " << ref_j);
	DEBUG(5, "free i " << pos_i << " free j " << pos_j);
	
	throw std::runtime_error("SHAKE failure in ??? (SHAKE)");
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
  skip_next.assign(skip_next.size(), true);
  
  return convergence;

}

