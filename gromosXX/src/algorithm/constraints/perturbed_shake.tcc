/**
 * @file perturbed_shake.tcc
 * contains the template methods for
 * the class Perturbed_Shake.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * Constructor.
 */
template<math::virial_enum do_virial>
algorithm::Perturbed_Shake<do_virial>
::Perturbed_Shake(Shake<do_virial> const & shake)
  : Algorithm("PerturbedShake"),
    m_shake(shake)
{
}

/**
 * Destructor.
 */
template<math::virial_enum do_virial>
algorithm::Perturbed_Shake<do_virial>
::~Perturbed_Shake()
{
}

/**
 * do one iteration
 */      
template<math::virial_enum do_virial, math::boundary_enum b>
static bool _perturbed_shake(topology::Topology const &topo,
		   configuration::Configuration & conf,
		   int const first,
		   std::vector<bool> &skip_now,
		   std::vector<bool> &skip_next,
		   std::vector<topology::perturbed_two_body_term_struct>
		   const & constr,
		   std::vector<interaction::bond_type_struct> const & param,
		   double const dt,
		   math::Periodicity<b> const & periodicity,
		   double const tolerance,
		   bool do_constraint_force = false,
		   size_t force_offset = 0)
{
  bool convergence = true;

  // index for constraint_force...
  size_t k = 0;
  double const dt2 = dt * dt;
  
  // and constraints
  for(typename std::vector<topology::perturbed_two_body_term_struct>
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
	
    double r0 = (1.0 - topo.lambda()) * param[it->A_type].r0 + 
      topo.lambda() * param[it->B_type].r0;
    
    double constr_length2 = r0 * r0;
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
      
      // the perturbed energy derivatives
      conf.current().perturbed_energy_derivatives.
	constraints_energy[topo.atom_energy_group()[it->i]] -=
	sqrt(dot(ref_r * lambda / dt2, ref_r * lambda / dt2)) *
	(param[it->B_type].r0 - param[it->A_type].r0);
      
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
 * shake perturbed solute
 */
template<math::virial_enum do_virial, math::boundary_enum b>
static int _perturbed_solute(topology::Topology const & topo,
			     configuration::Configuration & conf,
			     std::vector<interaction::bond_type_struct> 
			     const & param,
			     double dt, int const max_iterations,
			     double const tolerance)
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...

  DEBUG(8, "\tshaking perturbed SOLUTE");
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

    convergence = _perturbed_shake<do_virial, b>
      (topo, conf, first, skip_now, skip_next,
       topo.perturbed_solute().distance_constraints(), 
       param, dt,
       periodicity, tolerance, true);

    if(++num_iterations > max_iterations){
      io::messages.add("Perturbed SHAKE error. too many iterations",
		       "Perturbed_Shake::solute",
		       io::message::critical);
      throw std::runtime_error("SHAKE failure in perturbed solute");
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

  return num_iterations;

} // solute

/**
 * apply the SHAKE algorithm
 */
template<math::virial_enum do_virial>
int algorithm::Perturbed_Shake<do_virial>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "applying SHAKE");
  bool do_vel = false;
  int iterations = 0;
  
  // check whether we shake
  if (topo.perturbed_solute().distance_constraints().size() && 
      sim.param().shake.ntc > 1){
    DEBUG(8, "\twe need to shake perturbed SOLUTE");
    do_vel = true;
    switch(conf.boundary_type){
      case math::vacuum:
	iterations = _perturbed_solute<do_virial, math::vacuum>
	  (topo, conf, m_shake.parameter(), sim.time_step_size(), 
	   m_shake.max_iterations(), m_shake.tolerance());
	break;
      case math::rectangular:
	iterations = _perturbed_solute<do_virial, math::rectangular>
	  (topo, conf, m_shake.parameter(), sim.time_step_size(), 
	   m_shake.max_iterations(), m_shake.tolerance());
	break;
      case math::triclinic:
	iterations = _perturbed_solute<do_virial, math::triclinic>
	  (topo, conf, m_shake.parameter(), sim.time_step_size(),
	   m_shake.max_iterations(), m_shake.tolerance());
	break;
      default:
	throw std::string("wrong boundary type");
    }
  }
  
  // shaken velocity
  conf.current().vel = (conf.current().pos - conf.old().pos) / 
    sim.time_step_size();

  // return success!
  return 0;
		   
}

template<math::virial_enum do_virial>
int algorithm::Perturbed_Shake<do_virial>
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim)
{

  if (sim.param().start.shake_pos){
    std::cout << "shaking perturbed initial positions\n";

    // old and current pos and vel are the same...
    // shake the current ones
    apply(topo, conf, sim);

    // restore the velocities
    conf.current().vel = conf.old().vel;
    
    // take a step back
    conf.old().pos = conf.current().pos;
    
    if (sim.param().start.shake_vel){
      std::cout << "shaking initial velocities\n";

      conf.current().pos = conf.old().pos - 
	sim.time_step_size() * conf.old().vel;
    
      // shake again
      apply(topo, conf, sim);
    
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
