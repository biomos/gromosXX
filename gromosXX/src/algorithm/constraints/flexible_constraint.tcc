/**
 * @file flexible_constraint.tcc
 * implements flexible constraint algorithm
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

#include <util/debug.h>

/**
 * Constructor.
 */
template<math::virial_enum do_virial>
algorithm::Flexible_Constraint<do_virial>
::Flexible_Constraint(double const tolerance, int const max_iterations)
  : Algorithm("FlexibleShake"),
    m_tolerance(tolerance),
    m_max_iterations(max_iterations)
{
}

/**
 * Destructor
 */
template<math::virial_enum do_virial>
algorithm::Flexible_Constraint<do_virial>
::~Flexible_Constraint()
{
}

template<math::virial_enum do_virial>
void algorithm::Flexible_Constraint<do_virial>
::tolerance(double const tol)
{
  m_tolerance = tol;
}

/**
 * do one iteration
 */      
template<math::virial_enum do_virial, math::boundary_enum b>
static int _flexible_shake(topology::Topology const &topo,
			   configuration::Configuration & conf,
			   bool & convergence,
			   std::vector<double> & flex_len,
			   std::vector<bool> & skip_now,
			   std::vector<bool> & skip_next,
			   double const dt,
			   math::Periodicity<b> const & periodicity,
			   double const tolerance,
			   bool do_constraint_force = false)
{
  convergence = true;

  // index for constraint_force...
  size_t k = 0;
  double const dt2 = dt * dt;
  
  conf.special().flexible_ekin.assign(conf.special().flexible_ekin.size(), 0.0);
  conf.special().flexible_epot.assign(conf.special().flexible_epot.size(), 0.0);

  // and constraints
  for(typename std::vector<topology::two_body_term_struct>
	::const_iterator
	it = topo.solute().distance_constraints().begin(),
	to = topo.solute().distance_constraints().end();
      it != to;
      ++it, ++k ){
	
    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j]) continue;

    DEBUG(10, "i: " << it->i << " j: " << it->j);

    // the position
    math::Vec &pos_i = conf.current().pos(it->i);
    math::Vec &pos_j = conf.current().pos(it->j);

    DEBUG(10, "\ni: " << pos_i << "\nj: " << pos_j);
	
    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);
    double dist2 = dot(r, r);
	
    double constr_length2 = flex_len[k] * flex_len[k];
    double diff = constr_length2 - dist2;

    DEBUG(15, "constr: " << constr_length2 << " dist2: " << dist2);
	  
    if(fabs(diff) >= constr_length2 * tolerance * 2.0){
      // we have to shake
      DEBUG(10, "flexible shaking");
      
      // the reference position
      const math::Vec &ref_i = conf.old().pos(it->i);
      const math::Vec &ref_j = conf.old().pos(it->j);
      
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

	std::cout << "FLEXIBLE SHAKE ERROR (orthogonal vectors)\n"
		  << "\tatom i    : " << it->i << "\n"
		  << "\tatom j    : " << it->j << "\n"
		  << "\tref i     : " << math::v2s(ref_i) << "\n"
		  << "\tref j     : " << math::v2s(ref_j) << "\n"
		  << "\tfree i    : " << math::v2s(pos_i) << "\n"
		  << "\tfree j    : " << math::v2s(pos_j) << "\n"
		  << "\tref r     : " << math::v2s(ref_r) << "\n"
		  << "\tr         : " << math::v2s(r) << "\n"
		  << "\tsp        : " << sp << "\n"
		  << "\tconstr    : " << constr_length2 << "\n"
		  << "\tdiff      : " << diff << "\n"
		  << "\tforce i   : " << math::v2s(conf.old().force(it->i)) << "\n"
		  << "\tforce j   : " << math::v2s(conf.old().force(it->j)) << "\n"
		  << "\tvel i     : " << math::v2s(conf.current().vel(it->i)) << "\n"
		  << "\tvel j     : " << math::v2s(conf.current().vel(it->j)) << "\n"
		  << "\told vel i : " << math::v2s(conf.old().vel(it->i)) << "\n"
		  << "\told vel j : " << math::v2s(conf.old().vel(it->j)) << "\n\n";
	
	return E_SHAKE_FAILURE;
      }
	  
      // lagrange multiplier
      double lambda = diff / (sp * 2 *
			      (1.0 / topo.mass()(it->i) +
			       1.0 / topo.mass()(it->j) ));      

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
      pos_i += ref_r / topo.mass()(it->i);
      pos_j -= ref_r / topo.mass()(it->j);
	  
      convergence = false;

      // consider atoms in the next step
      skip_next[it->i] = false;
      skip_next[it->j] = false;
      
    } // we have to shake
  } // constraints
      
  
  return 0;

}    

template<math::virial_enum do_virial, math::boundary_enum b>
static void _calc_distance(topology::Topology const &topo,
			   configuration::Configuration & conf,
			   simulation::Simulation const & sim,
			   std::vector<interaction::bond_type_struct> const & param,
			   std::vector<double> & flex_len, double const dt)
{
  flex_len.clear();

  math::Periodicity<b> periodicity(conf.current().box);

  //loop over all constraints
  size_t k = 0;
  size_t com, ir;
  
  for(std::vector<topology::two_body_term_struct>::const_iterator
	it = topo.solute().distance_constraints().begin(),
	to = topo.solute().distance_constraints().end();
      it != to;
      ++it, ++k){
    
    // the position
    assert(conf.current().pos.size() > int(it->i));
    assert(conf.current().pos.size() > int(it->j));
    
    math::Vec const & pos_i = conf.current().pos(it->i);
    math::Vec const & pos_j = conf.current().pos(it->j);
	
    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);

    // unconstrained distance
    const double dist2 = dot(r, r);
    
    assert(conf.old().pos.size() > int(it->i));
    assert(conf.old().pos.size() > int(it->j));

    const math::Vec &ref_i = conf.old().pos(it->i);
    const math::Vec &ref_j = conf.old().pos(it->j);
      
    math::Vec ref_r;
    periodicity.nearest_image(ref_i, ref_j, ref_r);

    // reference distance
    const double ref_dist2 = dot(ref_r, ref_r);

    // standard formula with velocity along contsraints correction
    // (not the velocityless formula):
    assert(topo.mass().size() > int(it->i));
    assert(topo.mass().size() > int(it->j));

    const double red_mass = 1 / (1/topo.mass()(it->i) + 1/topo.mass()(it->j));
    const double dt2 = dt * dt;
      
    // calculate the force on constraint k
    assert(conf.special().flexible_vel.size() > k);

    const double force_on_constraint  = (red_mass / dt2) * 
      (sqrt(dist2) - sqrt(ref_dist2) - conf.special().flexible_vel[k] * dt);

    // zero energy distance

    // =================================================
    // it->type:  -> ideal bond length
    // flex_len:  flexible constraint distance
    // =================================================

    // const double constr_length2 = param(it->type).r0 * param(it->type).r0;
    
    // calculate the flexible constraint distance
    const double new_len = force_on_constraint / param[it->type].K + 
      param[it->type].r0;
    
    // store for shake
    flex_len.push_back(new_len);

    // update the velocity array
    conf.special().flexible_vel[k] = (new_len - sqrt(ref_dist2)) / dt;

    // now we have to store the kinetic energy of the constraint
    // length change in the correct temperature bath...
    sim.multibath().in_bath(it->i, com, ir);

    conf.special().flexible_ekin[ir] +=
      0.5 * red_mass * conf.special().flexible_vel[k] *
      conf.special().flexible_vel[k];

    DEBUG(8, "constr " << it->i << " bath " << ir << " ekin " 
	  << 0.5 * red_mass * conf.special().flexible_vel[k] * 
	  conf.special().flexible_vel[k]);


    // calculate Epot in the bond length constraints
    conf.special().flexible_epot[topo.atom_energy_group()[it->i]] += 
      0.5 * param[it->type].K * (param[it->type].r0 - new_len) * 
      (param[it->type].r0 - new_len);
      
    DEBUG(5, "flex_constraint_distance: " << new_len);
    
  } // loop over all constraints
  
}

/**
 * shake solute
 */
template<math::virial_enum do_virial, math::boundary_enum b>
static int _flexible_solute(topology::Topology const & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    std::vector<interaction::bond_type_struct> const & param,
			    double dt, int const max_iterations,
			    double const tolerance,
			    double & timing)
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...

  DEBUG(8, "\tFlexible shaking SOLUTE");
  math::Periodicity<b> periodicity(conf.current().box);
  
  std::vector<double> flex_len;
  _calc_distance<do_virial, b>(topo, conf, sim, param, flex_len, dt);

  // conf.constraint_force() = 0.0;
  // m_lambda = 0.0;

  const double start = util::now();

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);
  
  bool convergence = false;
  while(!convergence){
    DEBUG(9, "\titeration" << std::setw(10) << num_iterations);

    if(_flexible_shake<do_virial, b>
       (topo, conf, convergence, flex_len, skip_now, skip_next,
	dt,
	periodicity, tolerance, true)
       ){
      io::messages.add("Flexible SHAKE error. vectors orthogonal",
		       "Flexible_Constraint::solute",
		       io::message::error);
      std::cout << "Flexible SHAKE failure in solute!" << std::endl;
      return E_SHAKE_FAILURE_SOLUTE;
    }
    
    if(++num_iterations > max_iterations){
      io::messages.add("Flexible SHAKE error. too many iterations",
		       "Flexible_Constraint::solute",
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

  timing += util::now() - start;

  return 0;
}


/**
 * apply the Flexible SHAKE algorithm
 */
template<math::virial_enum do_virial>
int algorithm::Flexible_Constraint<do_virial>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "applying Flexible SHAKE");
  bool do_vel = false;
  int error = 0;
  
  // check whether we shake
  if (topo.solute().distance_constraints().size() && 
      sim.param().constraint.solute.algorithm == simulation::constr_flexshake &&
      sim.param().constraint.ntc > 1){
    DEBUG(8, "\twe need to flexible shake SOLUTE");
    do_vel = true;
    switch(conf.boundary_type){
      case math::vacuum:
	error = _flexible_solute<do_virial, math::vacuum>
	  (topo, conf, sim, parameter(), sim.time_step_size(), 
	   m_max_iterations, m_tolerance, m_timing);
	break;
      case math::rectangular:
	error = _flexible_solute<do_virial, math::rectangular>
	  (topo, conf, sim, parameter(), sim.time_step_size(), 
	   m_max_iterations, m_tolerance, m_timing);
	break;
      case math::triclinic:
	error = _flexible_solute<do_virial, math::triclinic>
	  (topo, conf, sim, parameter(), sim.time_step_size(),
	   m_max_iterations, m_tolerance, m_timing);
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

  // shaken velocity
  conf.current().vel = (conf.current().pos - conf.old().pos) / 
    sim.time_step_size();
  
  // return success!
  return 0;
		   
}

template<math::virial_enum do_virial>
int algorithm::Flexible_Constraint<do_virial>
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim)
{

  if (sim.param().start.shake_pos){
    std::cout << "(flexible) shaking initial positions\n";

    // old and current pos and vel are the same...
    conf.old().pos = conf.current().pos;
    conf.old().vel = conf.current().vel;

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
