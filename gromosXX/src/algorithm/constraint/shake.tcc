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

/**
 * Destructor.
 */
template<typename t_simulation>
algorithm::Shake<t_simulation>
::~Shake()
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
  
  sys.constraint_force() = 0.0;
  m_lambda = 0.0;

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  int first = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);
  
  bool convergence = false;
  while(!convergence){

    convergence = _shake(topo, sys, first, skip_now, skip_next,
			 topo.solute().distance_constraints(), true);

    if(++num_iterations > max_iterations){
      io::messages.add("SHAKE error. too many iterations",
		       "Shake::solute",
		       io::message::critical);
      throw std::runtime_error("SHAKE failure in solute");
    }

    skip_now = skip_next;
    skip_next.assign(skip_next.size(), true);

  } // convergence?

  // print out the constraint_force
  for (unsigned int i=0; i < topo.solute().distance_constraints().size();++i){
    sys.constraint_force()(i) *= 1 /(dt * dt);
    DEBUG(5, "constraint_force " << sqrt(dot(sys.constraint_force()(i),sys.constraint_force()(i)) ));
  }
  
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
			     topo.solvent(i).distance_constraints(), false);
	
	// std::cout << num_iterations+1 << std::endl;
	if(++num_iterations > max_iterations){
	  io::messages.add("SHAKE error. too many iterations",
			   "Shake::solvent",
			   io::message::critical);
	  throw std::runtime_error("SHAKE failure in solvent");
	}

	skip_now = skip_next;
	skip_next.assign(skip_next.size(), true);

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
template<typename t_distance_struct>
bool algorithm::Shake<t_simulation>
::_shake(typename simulation_type::topology_type const &topo,
	 typename simulation_type::system_type &sys,
	 int const first,
	 std::vector<bool> &skip_now,
	 std::vector<bool> &skip_next,
	 std::vector<t_distance_struct>
	 & constr, bool do_constraint_force, size_t force_offset)
{
  bool convergence = true;

  // index for constraint_force...
  size_t k = 0;
  
  // and constraints
  for(typename std::vector<t_distance_struct>
	::const_iterator
	it = constr.begin(),
	to = constr.end();
      it != to;
      ++it, ++k ){
	
    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j]) continue;

    DEBUG(10, "i: " << it->i << " j: " << it->j << " first: " << first);

    // the position
    math::Vec &pos_i = sys.pos()(first+it->i);
    math::Vec &pos_j = sys.pos()(first+it->j);

    DEBUG(10, "\ni: " << pos_i << "\nj: " << pos_j);
	
    math::Vec r;
    sys.periodicity().nearest_image(pos_i, pos_j, r);
    double dist2 = dot(r, r);
	
    double constr_length2 = it->b0 * it->b0;
    double diff = constr_length2 - dist2;

    DEBUG(15, "constr: " << constr_length2 << " dist2: " << dist2);
	  
    if(fabs(diff) >= constr_length2 * m_tolerance * 2.0){
      // we have to shake
      DEBUG(10, "shaking");
      
      // the reference position
      const math::Vec &ref_i = sys.old_pos()(first+it->i);
      const math::Vec &ref_j = sys.old_pos()(first+it->j);
      
      math::Vec ref_r;
      sys.periodicity().nearest_image(ref_i, ref_j, ref_r);

      double sp = dot(ref_r, r);
	  
      if(sp < constr_length2 * math::epsilon){
	io::messages.add("SHAKE error. vectors orthogonal",
			 "Shake::???",
			 io::message::critical);
	DEBUG(5, "ref i " << ref_i << " ref j " << ref_j);
	DEBUG(5, "free i " << pos_i << " free j " << pos_j);
	DEBUG(5, "ref r " << ref_r);
	DEBUG(5, "r " << r);
	
	throw std::runtime_error("SHAKE failure in ??? (SHAKE)");
      }
	  
      // lagrange multiplier
      double lambda = diff / (sp * 2 *
			      (1.0 / topo.mass()(first+it->i) +
			       1.0 / topo.mass()(first+it->j) ));      

      DEBUG(10, "lagrange multiplier " << lambda);

      if (do_constraint_force == true){
	
	//if it is a solute sum up constraint forces
	assert(unsigned(sys.constraint_force().size()) > k + force_offset);
	sys.constraint_force()(k+force_offset) += (lambda * ref_r);
	m_lambda(k) += lambda;
      }

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
      
  
  return convergence;

}

/**
 * add bond type.
 */
template<typename t_simulation>
inline void algorithm::Shake<t_simulation>
::add_bond_type(interaction::bond_type_struct s)
{
  m_bond_parameter.push_back(s);
}

/**
 * add bond type.
 */
template<typename t_simulation>
inline void algorithm::Shake<t_simulation>
::add_bond_type(double K, double r0)
{
  interaction::bond_type_struct s;
  s.K = K;
  s.r0 = r0;
  add_bond_type(s);
}

/**
 * add all bonds to the solute constraint vector and
 * remove them from the bond vector.
 */
template<typename t_simulation>
inline void 
algorithm::Shake<t_simulation>
::add_bond_length_constraints(typename t_simulation::topology_type &topo)
{
  simulation::Solute & solute = topo.solute();

  std::vector<simulation::Bond> bonds;
  std::vector<simulation::Bond>::iterator it = solute.bonds().begin(),
    to = solute.bonds().end();
  for( ; it != to; ++it){
    DEBUG(9, "constraint " << it->i << " - " << it->j);
    solute.add_distance_constraint(it->i, it->j, m_bond_parameter[it->type].r0);
  }
  solute.bonds() = bonds;
}

/**
 * add bonds connecting an atom of type iac to the
 * constraint vector and remove from the bond vector.
 */
template<typename t_simulation>
inline void
algorithm::Shake<t_simulation>
::add_bond_length_constraints(int iac,
			      std::vector<int> const &atom_iac,
			      typename t_simulation::topology_type &topo)
{
  simulation::Solute & solute = topo.solute();

  std::vector<simulation::Bond> bonds;
  std::vector<simulation::Bond>::iterator it = solute.bonds().begin(),
    to = solute.bonds().end();

  for( ; it != to; ++it){
    if(atom_iac[it->i] == iac || atom_iac[it->j] == iac){
      DEBUG(9, "constraint " << it->i << " - " << it->j);      
      solute.add_distance_constraint(it->i, it->j, m_bond_parameter[it->type].r0);
    }
    else
      bonds.push_back(simulation::Bond(it->i, it->j, it->type));
  }
  solute.bonds() = bonds;
}
    
/**
 * add bonds connecting an atom of mass mass to the
 * constraint vector and remove from the bond vector.
 */
template<typename t_simulation>
inline void
algorithm::Shake<t_simulation>
::add_bond_length_constraints(double mass,
			      math::SArray const &atom_mass,
			      typename t_simulation::topology_type &topo)
{
  simulation::Solute & solute = topo.solute();
  
  std::vector<simulation::Bond> bonds;
  std::vector<simulation::Bond>::iterator it = solute.bonds().begin(),
    to = solute.bonds().end();

  for( ; it != to; ++it){
    if(atom_mass(it->i) == mass || atom_mass(it->j) == mass){
      DEBUG(9, "constraint " << it->i << " - " << it->j);      
      solute.add_distance_constraint(it->i, it->j, m_bond_parameter[it->type].r0);
    }
    else
      bonds.push_back(simulation::Bond(it->i, it->j, it->type));
  }
  solute.bonds() = bonds;
}

template<typename t_simulation>
inline void
algorithm::Shake<t_simulation>
::init(t_simulation &sim, io::Argument &args, io::InTopology &topo,
       io::InInput &input)
{
  // initialize SHAKE / ??
  DEBUG(7, "md init shake");
  int ntc;
  double tolerance;
  input.read_SHAKE(ntc, tolerance);
  this->tolerance(tolerance);

  std::cout << "DISTANCE CONSTRAINTS\n";

  switch(ntc){
    case 1:
      break;
    case 2: 
      std::cout << "\tSHAKE bonds containing hydrogen atoms" << std::endl;
      // read in the parameter
      topo >> *this;
      add_bond_length_constraints(1.0,
				  sim.topology().mass(),
				  sim.topology());
      break;
    case 3: 
      std::cout << "\tSHAKE all bonds" << std::endl;
      // read in parameter
      topo >> *this;
      add_bond_length_constraints(sim.topology());
      break;
    default:
      std::cout << "\twrong ntc parameter" << std::endl;
  }

  // give the constraint force the correct size...
  sim.system().constraint_force().resize(sim.topology().solute().
				distance_constraints().size());
  m_lambda.resize(sim.topology().solute().distance_constraints().size());

  std::cout << "END\n";

}


