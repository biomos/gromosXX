/**
 * @file shake.tcc
 * contains the template methods for
 * the class Shake.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

#include <util/error.h>

/**
 * Constructor.
 */
template<math::virial_enum do_virial>
algorithm::Shake<do_virial>
::Shake(double const tolerance, int const max_iterations)
  : Algorithm("Shake"),
    m_tolerance(tolerance),
    m_max_iterations(max_iterations)
{
}

/**
 * Destructor.
 */
template<math::virial_enum do_virial>
algorithm::Shake<do_virial>
::~Shake()
{
}

template<math::virial_enum do_virial>
void algorithm::Shake<do_virial>
::tolerance(double const tol)
{
  m_tolerance = tol;
}

/**
 * do one iteration
 */      
template<math::virial_enum do_virial, math::boundary_enum b>
static int _shake(topology::Topology const &topo,
		  configuration::Configuration & conf,
		  bool & convergence,
		  int const first,
		  std::vector<bool> &skip_now,
		  std::vector<bool> &skip_next,
		  std::vector<topology::two_body_term_struct> const & constr,
		  std::vector<interaction::bond_type_struct> const & param,
		  double const dt,
		  math::Periodicity<b> const & periodicity,
		  double const tolerance,
		  bool do_constraint_force = false, size_t force_offset = 0)
{
  convergence = true;

  // index for constraint_force...
  size_t k = 0;
  double const dt2 = dt * dt;
  
  // and constraints
  for(typename std::vector<topology::two_body_term_struct>
	::const_iterator
	it = constr.begin(),
	to = constr.end();
      it != to;
      ++it, ++k ){
	
    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j]) continue;

    DEBUG(10, "i: " << it->i << " j: " << it->j << " first: " << first);

    // the position
    math::Vec &pos_i = conf.current().pos(first+it->i);
    math::Vec &pos_j = conf.current().pos(first+it->j);

    DEBUG(10, "\ni: " << pos_i << "\nj: " << pos_j);
	
    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);
    double dist2 = dot(r, r);
	
    double constr_length2 = param[it->type].r0 * param[it->type].r0;
    double diff = constr_length2 - dist2;

    DEBUG(15, "constr: " << constr_length2 << " dist2: " << dist2);
	  
    if(fabs(diff) >= constr_length2 * tolerance * 2.0){
      // we have to shake
      DEBUG(10, "shaking");
      
      // the reference position
      const math::Vec &ref_i = conf.old().pos(first+it->i);
      const math::Vec &ref_j = conf.old().pos(first+it->j);
      
      math::Vec ref_r;
      periodicity.nearest_image(ref_i, ref_j, ref_r);

      double sp = dot(ref_r, r);
	  
      if(sp < constr_length2 * math::epsilon){
	/*
	io::messages.add("SHAKE error. vectors orthogonal",
			 "Shake::???",
			 io::message::critical);
	*/
	DEBUG(5, "ref i " << ref_i << " ref j " << ref_j);
	DEBUG(5, "free i " << pos_i << " free j " << pos_j);
	DEBUG(5, "ref r " << ref_r);
	DEBUG(5, "r " << r);

	std::cout << "SHAKE ERROR\n"
		  << "\tatom i    : " << it->i << "\n"
		  << "\tatom j    : " << it->j << "\n"
		  << "\tfirst     : " << first << "\n"
		  << "\tref i     : " << math::v2s(ref_i) << "\n"
		  << "\tref j     : " << math::v2s(ref_j) << "\n"
		  << "\tfree i    : " << math::v2s(pos_i) << "\n"
		  << "\tfree j    : " << math::v2s(pos_j) << "\n"
		  << "\tref r     : " << math::v2s(ref_r) << "\n"
		  << "\tr         : " << math::v2s(r) << "\n"
		  << "\tsp        : " << sp << "\n"
		  << "\tconstr    : " << constr_length2 << "\n"
		  << "\tdiff      : " << diff << "\n"
		  << "\tforce i   : " << math::v2s(conf.old().force(first+it->i)) << "\n"
		  << "\tforce j   : " << math::v2s(conf.old().force(first+it->j)) << "\n"
		  << "\tvel i     : " << math::v2s(conf.current().vel(first+it->i)) << "\n"
		  << "\tvel j     : " << math::v2s(conf.current().vel(first+it->j)) << "\n"
		  << "\told vel i : " << math::v2s(conf.old().vel(first+it->i)) << "\n"
		  << "\told vel j : " << math::v2s(conf.old().vel(first+it->j)) << "\n\n";
	
	return E_SHAKE_FAILURE;
      }
	  
      // lagrange multiplier
      double lambda = diff / (sp * 2 *
			      (1.0 / topo.mass()(first+it->i) +
			       1.0 / topo.mass()(first+it->j) ));      

      DEBUG(10, "lagrange multiplier " << lambda);

      /*
      if (do_constraint_force == true){
	
	//if it is a solute sum up constraint forces
	assert(unsigned(sys.constraint_force().size()) > k + force_offset);
	sys.constraint_force()(k+force_offset) += (lambda * ref_r);
	// m_lambda(k) += lambda;
      }
      */

      if (do_virial == math::atomic_virial){
	for(int a=0; a<3; ++a){
	  for(int aa=0; aa<3; ++aa){
	    conf.current().virial_tensor(a,aa) +=
	      ref_r(a) * ref_r(aa) * lambda / dt2;
	  }
	}
	DEBUG(12, "\tatomic virial done");
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
      
  
  return 0;

}    

/**
 * shake solute
 */
template<math::virial_enum do_virial, math::boundary_enum b>
static int solute(topology::Topology const & topo,
		  configuration::Configuration & conf,
		  std::vector<interaction::bond_type_struct> const & param,
		  double dt, int const max_iterations,
		  double const tolerance)
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...

  DEBUG(8, "\tshaking SOLUTE");
  math::Periodicity<b> periodicity(conf.current().box);
  
  // conf.constraint_force() = 0.0;
  // m_lambda = 0.0;

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  int first = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);
  
  bool convergence = false;
  while(!convergence){
    DEBUG(9, "\titeration" << std::setw(10) << num_iterations);

    if(_shake<do_virial, b>
       (topo, conf, convergence, first, skip_now, skip_next,
	topo.solute().distance_constraints(), param, dt,
	periodicity, tolerance, true)
       ){
      io::messages.add("SHAKE error. vectors orthogonal",
		       "Shake::solute",
		       io::message::error);
      std::cout << "SHAKE failure in solute!" << std::endl;
      return E_SHAKE_FAILURE_SOLUTE;
    }
    
    if(++num_iterations > max_iterations){
      io::messages.add("SHAKE error. too many iterations",
		       "Shake::solute",
		       io::message::error);
      return E_SHAKE_FAILURE_SOLUTE;
    }

    skip_now = skip_next;
    skip_next.assign(skip_next.size(), true);

  } // convergence?

  // constraint_force
  /*
  for (unsigned int i=0; i < topo.solute().distance_constraints().size();++i){
    conf.constraint_force()(i) *= 1 /(dt * dt);
    DEBUG(5, "constraint_force " << sqrt(dot(conf.constraint_force()(i),
					     conf.constraint_force()(i)) ));
  }
  */

  return 0;

} // solute


/**
 * shake solvent.
 */
template<math::virial_enum do_virial, math::boundary_enum b>
static int solvent(topology::Topology const & topo,
		   configuration::Configuration & conf,
		   std::vector<interaction::bond_type_struct> &param,
		   double dt, int const max_iterations, 
		   double const tolerance)
{

  DEBUG(8, "\tshaking SOLVENT");
  
  // the first atom of a solvent
  size_t first = topo.num_solute_atoms();

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int tot_iterations = 0;

  math::Periodicity<b> periodicity(conf.current().box);

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
	DEBUG(9, "\titeration" << std::setw(10) << num_iterations);

	if(_shake<do_virial, b>
	   (topo, conf, convergence, first, skip_now, skip_next,
	    topo.solvent(i).distance_constraints(), param, dt,
	    periodicity, tolerance, false)){
	  
	  io::messages.add("SHAKE error. vectors orthogonal",
			   "Shake::solute", io::message::error);
	  
	  std::cout << "SHAKE failure in solute!" << std::endl;
	  return E_SHAKE_FAILURE_SOLVENT;
	}
	
	// std::cout << num_iterations+1 << std::endl;
	if(++num_iterations > max_iterations){
	  io::messages.add("SHAKE error. too many iterations",
			   "Shake::solvent",
			   io::message::critical);
	  // throw std::runtime_error("SHAKE failure in solvent");
	  return E_SHAKE_FAILURE_SOLVENT;
	}

	skip_now = skip_next;
	skip_next.assign(skip_next.size(), true);

      } // while(!convergence)
      
      tot_iterations += num_iterations;
      
    } // molecules
    
  } // solvents

  return 0;
  
} // shake solvent




/**
 * apply the SHAKE algorithm
 */
template<math::virial_enum do_virial>
int algorithm::Shake<do_virial>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "applying SHAKE");
  bool do_vel = false;
  int error = 0;
  
  // check whether we shake
  if (topo.solute().distance_constraints().size() && 
      sim.param().shake.ntc > 1){
    DEBUG(8, "\twe need to shake SOLUTE");
    do_vel = true;
    switch(conf.boundary_type){
      case math::vacuum:
	error = solute<do_virial, math::vacuum>
	  (topo, conf, parameter(), sim.time_step_size(), 
	   m_max_iterations, m_tolerance);
	break;
      case math::rectangular:
	error = solute<do_virial, math::rectangular>
	  (topo, conf, parameter(), sim.time_step_size(), 
	   m_max_iterations, m_tolerance);
	break;
      case math::triclinic:
	error = solute<do_virial, math::triclinic>
	  (topo, conf, parameter(), sim.time_step_size(),
	   m_max_iterations, m_tolerance);
	break;
      default:
	throw std::string("wrong boundary type");
	
    }
  }
  
  if (error){
    std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLUTE" << std::endl;
    // save old positions to final configuration... (even before free-flight!)
    conf.current().pos = conf.old().pos;
    return E_SHAKE_FAILURE_SOLUTE;
  }
  
  if (sim.param().system.nsm){
    DEBUG(8, "\twe need to shake SOLVENT");
    do_vel = true;
    switch(conf.boundary_type){
      case math::vacuum:
	error = 
	  solvent<do_virial, math::vacuum>
	  (topo, conf, parameter(), sim.time_step_size(),
	   m_max_iterations, m_tolerance);
	break;
      case math::rectangular:
	error = 
	  solvent<do_virial, math::rectangular>
	  (topo, conf, parameter(), sim.time_step_size(),
	   m_max_iterations, m_tolerance);
	break;
      case math::triclinic:
	error = 
	  solvent<do_virial, math::triclinic>
	  (topo, conf, parameter(), sim.time_step_size(),
	   m_max_iterations, m_tolerance);
	break;
      default:
	throw std::string("wrong boundary type");
	
    }
  }
  
  if (error){
    std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLVENT" << std::endl;
    // save old positions to final configuration... (even before free-flight!)
    conf.current().pos = conf.old().pos;
    return E_SHAKE_FAILURE_SOLVENT;
  }
      
  // shaken velocity
  conf.current().vel = (conf.current().pos - conf.old().pos) / 
    sim.time_step_size();

  // return success!
  return 0;
		   
}

template<math::virial_enum do_virial>
int algorithm::Shake<do_virial>
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim)
{

  if (sim.param().start.shake_pos){
    std::cout << "shaking initial positions\n";

    // old and current pos and vel are the same...
    // shake the current ones
    if (apply(topo, conf, sim))
      return E_SHAKE_FAILURE;

    // restore the velocities
    conf.current().vel = conf.old().vel;
    
    // take a step back
    conf.old().pos = conf.current().pos;
    
    if (sim.param().start.shake_vel){
      std::cout << "shaking initial velocities\n";

      conf.current().pos = conf.old().pos - 
	sim.time_step_size() * conf.old().vel;
    
      // shake again
      if (apply(topo, conf, sim))
	return E_SHAKE_FAILURE;
    
      // restore the positions
      conf.current().pos = conf.old().pos;
    
      // velocities are in opposite direction (in time)
      conf.current().vel = -1.0 * conf.current().vel;
      conf.old().vel = conf.current().vel;
    }
    
  }
  else if (sim.param().start.shake_vel){
    io::messages.add("shaking velocities without shaking positions illegal.",
		     "shake", io::message::error);
  }
  
  return 0;
}
