/**
 * @file flexible_constraint.cc
 * implements flexible constraint algorithm
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction.h>
#include <interaction/nonbonded/interaction/nonbonded_interaction.h>

#include <interaction/forcefield/forcefield.h>

#include <math/periodicity.h>

#include <algorithm/constraints/flexible_constraint.h>

#include <util/template_split.h>
#include <util/error.h>
#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * Constructor.
 */
algorithm::Flexible_Constraint
::Flexible_Constraint(double const tolerance, int const max_iterations,
		      interaction::Forcefield * ff)
  : Algorithm("FlexibleShake"),
    m_tolerance(tolerance),
    m_max_iterations(max_iterations)
{

  if (ff){
    for(size_t i=0; i<ff->size(); ++i){
      if ((*ff)[i]->name == "NonBonded"){
	// we have a nonbonded, try to cast it
	m_nonbonded = dynamic_cast<interaction::Nonbonded_Interaction *>((*ff)[i]);
      }
    }
    
    if(!ff){
      io::messages.add("Accessing the Nonbonded_Interaction failed!",
		       "Flexible_Constraint::Constructor",
		       io::message::error);
    }
  }
}

/**
 * Destructor
 */
algorithm::Flexible_Constraint
::~Flexible_Constraint()
{
}

void algorithm::Flexible_Constraint
::tolerance(double const tol)
{
  m_tolerance = tol;
}

/**
 * apply the Flexible SHAKE algorithm
 */
int algorithm::Flexible_Constraint
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "applying Flexible SHAKE");
  bool do_vel = false;
  int error = 0;

  conf.special().flexible_constraint.flexible_ekin.assign
    (conf.special().flexible_constraint.flexible_ekin.size(), 0.0);

  // check whether we shake
  if (topo.solute().distance_constraints().size() && 
      sim.param().constraint.solute.algorithm == simulation::constr_flexshake &&
      sim.param().constraint.ntc > 1){
    DEBUG(8, "\twe need to flexible shake SOLUTE");
    do_vel = true;

    SPLIT_VIRIAL_BOUNDARY(solute,
			  topo, conf, sim, m_force, error);
  }
  
  if (error){
    std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLUTE" << std::endl;
    // save old positions to final configuration... (even before free-flight!)
    conf.current().pos = conf.old().pos;
    return E_SHAKE_FAILURE_SOLUTE;
  }

  // return success!
  return 0;
}


/**
 * do one iteration
 */      
template<math::boundary_enum B, math::virial_enum V>
int algorithm::Flexible_Constraint::iteration
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 bool & convergence,
 std::vector<double> & flex_len,
 std::vector<bool> & skip_now,
 std::vector<bool> & skip_next,
 double dt,
 math::Periodicity<B> const & periodicity,
 bool do_constraint_force
 )
{
  convergence = true;

  // index for constraint_force...
  unsigned int k = 0;
  double const dt2 = dt * dt;
  
  // and constraints
  for(typename std::vector<topology::two_body_term_struct>::const_iterator
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

    DEBUG(10, "\ni: " << math::v2s(pos_i) << "\nj: " << math::v2s(pos_j));
	
    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);
    double dist2 = math::abs2(r);
	
    double constr_length2 = flex_len[k] * flex_len[k];
    double diff = constr_length2 - dist2;

    DEBUG(15, "constr: " << constr_length2 << " dist2: " << dist2);
	  
    if(fabs(diff) >= constr_length2 * m_tolerance * 2.0){
      // we have to shake
      DEBUG(10, "flexible shaking");
      
      // the reference position
      const math::Vec &ref_i = conf.old().pos(it->i);
      const math::Vec &ref_j = conf.old().pos(it->j);
      
      math::Vec ref_r;
      periodicity.nearest_image(ref_i, ref_j, ref_r);

      double sp = math::dot(ref_r, r);
	  
      if(sp < constr_length2 * math::epsilon){
	DEBUG(5, "ref i " << math::v2s(ref_i) << " ref j " << math::v2s(ref_j));
	DEBUG(5, "free i " << math::v2s(pos_i) << " free j " << math::v2s(pos_j));
	DEBUG(5, "ref r " << math::v2s(ref_r));
	DEBUG(5, "r " << math::v2s(r));

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

      if (V == math::atomic_virial){
	for(int a=0; a<3; ++a){
	  for(int aa=0; aa<3; ++aa){
	    conf.old().virial_tensor(a,aa) +=
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

//================================================================================
//================================================================================
// 
// exact interation
//
//================================================================================
//================================================================================

template<math::boundary_enum B, math::virial_enum V>
int algorithm::Flexible_Constraint::exact_iteration
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 bool & convergence,
 std::vector<double> & flex_len,
 std::vector<std::vector<math::Vec> > & force,
 double dt,
 math::Periodicity<B> const & periodicity
 )
{

  convergence = true;

  unsigned int k = 0;
  double const dt2 = dt * dt;
  
  // and constraints
  for(typename std::vector<topology::two_body_term_struct>::const_iterator
	it = topo.solute().distance_constraints().begin(),
	to = topo.solute().distance_constraints().end();
      it != to;
      ++it, ++k ){
	
    // no skipping for the exact algorithm
    // ('cause every constraint changes ALL positions)

    DEBUG(8, "i: " << it->i << " j: " << it->j);

    // the position
    math::Vec &pos_i = conf.current().pos(it->i);
    math::Vec &pos_j = conf.current().pos(it->j);

    DEBUG(10, "\ni: " << math::v2s(pos_i) << "\nj: " << math::v2s(pos_j));
    
    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);
    double dist2 = math::abs2(r);
    
    double constr_length2 = flex_len[k] * flex_len[k];
    double diff = constr_length2 - dist2;

    DEBUG(12, "constr: " << constr_length2 << " dist2: " << dist2);
    
    if(fabs(diff) >= constr_length2 * m_tolerance * 2.0){
      // we have to shake
      DEBUG(10, "exact flexible shaking");
      
      // the undetermined constraint forces
      // (as function of the reference position...)
      const math::Vec &f_i = force[k][it->i];
      const math::Vec &f_j = force[k][it->j];
      
      math::Vec ud_f = 
	1.0 / topo.mass()(it->i) * f_i -
	1.0 / topo.mass()(it->j) * f_j;

      double sp = dot(r, ud_f);
	  
      if(sp < constr_length2 * math::epsilon){

	DEBUG(5, "f(ud) i " << math::v2s(f_i) << " f(ud) j " << math::v2s(f_j));
	DEBUG(5, "free i " << math::v2s(pos_i) << " free j " << math::v2s(pos_j));
	DEBUG(5, "ud_f " << math::v2s(ud_f));
	DEBUG(5, "r " << math::v2s(r));

	std::cout << "EXACT FLEXIBLE SHAKE ERROR (orthogonal vectors)\n"
		  << "\tatom i    : " << it->i << "\n"
		  << "\tatom j    : " << it->j << "\n"
		  << "\tf i       : " << math::v2s(f_i) << "\n"
		  << "\tf j       : " << math::v2s(f_j) << "\n"
		  << "\tfree i    : " << math::v2s(pos_i) << "\n"
		  << "\tfree j    : " << math::v2s(pos_j) << "\n"
		  << "\tud f      : " << math::v2s(ud_f) << "\n"
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
      const double lambda = diff / (- 4 * dt2 * sp);
      DEBUG(10, "lagrange multiplier " << lambda);

      if (V == math::atomic_virial){
	/*
	for(int a=0; a<3; ++a){
	  for(int aa=0; aa<3; ++aa){
	    conf.old().virial_tensor(a,aa) +=
	      ref_r(a) * ref_r(aa) * lambda / dt2;
	  }
	}
	DEBUG(12, "\tatomic virial done");
	*/
	// not implemented
	assert(false);
      }
      
      // update positions
      for(unsigned int a=0; a<topo.num_atoms(); ++a){
	conf.current().pos(a) -= 2 * dt2 / topo.mass()(a) * lambda * force[k][a];
      }
	  
      convergence = false;
      
    } // we have to shake
  } // constraints
      
  
  return 0;

}    


//================================================================================
// flexible constraint distance
//================================================================================

template<math::boundary_enum B, math::virial_enum V>
void algorithm::Flexible_Constraint::calc_distance
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim,
 std::vector<double> & flex_len
 )
{

  DEBUG(7, "FLEXIBLE CONSTRAINTS::calc_distance");
  
  flex_len.clear();

  math::Periodicity<B> periodicity(conf.current().box);

  //loop over all constraints
  unsigned int k = 0;
  unsigned int com, ir;
  
  for(std::vector<topology::two_body_term_struct>::const_iterator
	it = topo.solute().distance_constraints().begin(),
	to = topo.solute().distance_constraints().end();
      it != to;
      ++it, ++k){
    
    DEBUG(8, "flexible constraint " << k);

    assert(conf.old().pos.size() > (it->i));
    assert(conf.old().pos.size() > (it->j));

    const math::Vec &ref_i = conf.old().pos(it->i);
    const math::Vec &ref_j = conf.old().pos(it->j);
      
    math::Vec ref_r;
    periodicity.nearest_image(ref_i, ref_j, ref_r);

    // reference distance
    const double ref_dist2 = math::abs2(ref_r);

    assert(topo.mass().size() > (it->i));
    assert(topo.mass().size() > (it->j));

    const double red_mass = 1 / (1/topo.mass()(it->i) + 1/topo.mass()(it->j));
    const double dt2 = sim.time_step_size() * sim.time_step_size();

    DEBUG(10, "red mass = " << red_mass << "   dt2 = " << dt2);
      
    // calculate the force on constraint k
    double force_on_constraint;

    //**************************************************
    // velocities along constraints ???
    //**************************************************
    
    if (sim.param().constraint.solute.flexshake_mode == 0 ||
	sim.param().constraint.solute.flexshake_mode == 2){
      // the position
      assert(conf.current().pos.size() > (it->i));
      assert(conf.current().pos.size() > (it->j));
      
      math::Vec const & pos_i = conf.current().pos(it->i);
      math::Vec const & pos_j = conf.current().pos(it->j);
      
      math::Vec r;
      periodicity.nearest_image(pos_i, pos_j, r);
      
      // unconstrained distance
      const double dist2 = math::abs2(r);
      DEBUG(10, "unconstrained distance = " << dist2);
      
      // take out velocities along constraint from previous step...
      assert(conf.special().flexible_constraint.flexible_vel.size() > k);
      DEBUG(10, "flexible constraint velocity =  "
	    << conf.special().flexible_constraint.flexible_vel[k]);

      force_on_constraint = (red_mass / dt2) * 
	(sqrt(dist2) - sqrt(ref_dist2) - 
	 conf.special().flexible_constraint.flexible_vel[k] * sim.time_step_size());
    }
    else{
      // only the potenital energy
      // or in other words
      // the forces on the constraint, no velocities
      math::Vec pos_i = conf.old().pos(it->i) +
	dt2 / topo.mass()(it->i) * conf.old().force(it->i);
      math::Vec pos_j = conf.old().pos(it->j) +
	dt2 / topo.mass()(it->j) * conf.old().force(it->j);

      math::Vec r;
      periodicity.nearest_image(pos_i, pos_j, r);

      // unconstrained distance
      const double dist2 = math::abs2(r);
      DEBUG(10, "unconstrained distance = " << dist2);

      // ignore velocities along constraints
      force_on_constraint = (red_mass / dt2) * 
	(sqrt(dist2) - sqrt(ref_dist2));
    }
      
    // zero energy distance

    // =================================================
    // it->type:  -> ideal bond length
    // flex_len:  flexible constraint distance
    // =================================================

    // const double constr_length2 = m_parameter(it->type).r0 * m_parameter(it->type).r0;
    
    // calculate the flexible constraint distance
    DEBUG(10, "F(c) = " << force_on_constraint);
    DEBUG(10, "K = " << m_parameter[it->type].K << "     r0 = " << m_parameter[it->type].r0);

    const double new_len = force_on_constraint / m_parameter[it->type].K + 
      m_parameter[it->type].r0;
    
    // store for shake
    flex_len.push_back(new_len);

    // update the velocity array
    conf.special().flexible_constraint.flexible_vel[k] = (new_len - sqrt(ref_dist2)) / 
      sim.time_step_size();

    // now we have to store the kinetic energy of the constraint
    // length change in the correct temperature bath...
    sim.multibath().in_bath(it->i, com, ir);

    conf.special().flexible_constraint.flexible_ekin[ir] +=
      0.5 * red_mass * conf.special().flexible_constraint.flexible_vel[k] *
      conf.special().flexible_constraint.flexible_vel[k];

    DEBUG(8, "constr " << it->i << " bath " << ir << " ekin " 
	  << 0.5 * red_mass * conf.special().flexible_constraint.flexible_vel[k] * 
	  conf.special().flexible_constraint.flexible_vel[k]);

    // calculate Epot in the bond length constraints
    conf.old().energies.constraints_energy[topo.atom_energy_group()[it->i]] += 
      0.5 * m_parameter[it->type].K * (m_parameter[it->type].r0 - new_len) * 
      (m_parameter[it->type].r0 - new_len);
      
    DEBUG(5, "flex_constraint_distance: " << new_len);
    
  } // loop over all constraints
}

//================================================================================
// undetermined forces
//================================================================================

template<math::boundary_enum B, math::virial_enum V>
void algorithm::Flexible_Constraint::calc_undetermined_forces
(
 topology::Topology &topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::vector<double> & flex_len, 
 std::vector<std::vector<math::Vec> > & force
 )
{
  DEBUG(7, "FLEXIBLE CONSTRAINTS: _calc_undetermined_forces");
  
  // resize to number of constraints
  force.resize(topo.solute().distance_constraints().size());

  math::Periodicity<B> periodicity(conf.current().box);

  const double dt2 = sim.time_step_size() * sim.time_step_size();

  //loop over all constraints
  unsigned int k = 0;
  for(std::vector<topology::two_body_term_struct>::const_iterator
	it = topo.solute().distance_constraints().begin(),
	to = topo.solute().distance_constraints().end();
      it != to;
      ++it, ++k){

    DEBUG(8, "\tconstraint " << k);
    
    force[k].assign(topo.num_atoms(), math::Vec(0.0));

    for(unsigned int a=0; a<topo.num_atoms(); ++a){

      assert(conf.current().pos.size() > (it->i));
      assert(conf.current().pos.size() > (it->j));
    
      // r_k^nc(t+dt) vector
      math::Vec r_nc;

      if (sim.param().constraint.solute.flexshake_mode == 2){
	// with velocities => use conf.current()
        math::Vec const & pos_i = conf.current().pos(it->i);
        math::Vec const & pos_j = conf.current().pos(it->j);
      
        periodicity.nearest_image(pos_i, pos_j, r_nc);
      }
      else{
	// no velocities => recalculate current positions from forces only
	math::Vec pos_i = conf.old().pos(it->i) +
	  dt2 * conf.old().force(it->i) / topo.mass()(it->i);
	math::Vec pos_j = conf.old().pos(it->j) +
	  dt2 * conf.old().force(it->j) / topo.mass()(it->j);
	
	periodicity.nearest_image(pos_i, pos_j, r_nc);
      }

      DEBUG(10, "r_nc = " << math::v2s(r_nc));

      // unconstrained distance
      const double dist = sqrt(math::abs2(r_nc));
      // normalise
      r_nc /= dist;
      
      // flexible constraint length
      assert(flex_len.size() > k);
      const double dk = flex_len[k];
      // reduced mass
      const double mu = 1.0 / (1.0 / topo.mass()(it->i) + 1.0 / topo.mass()(it->j));

      if (a == it->i){
	// constraint forces on atom i (k1) of constraint k : Eq 26
	DEBUG(12, "\tforces on constraint - atoms");

	math::Vec const & ref_i = conf.old().pos(it->i);
	math::Vec const & ref_j = conf.old().pos(it->j);
      
	// r_k(t) vector
	math::Vec r;
	periodicity.nearest_image(ref_i, ref_j, r);

	// last step (constrained) distance
	const double r_dist = sqrt(math::abs2(r));

	math::Matrix h1(0.0), h2(0.0), u1(0.0), u2(0.0);
	
	// loop over all other atoms
	for(unsigned int a2=0; a2<topo.num_atoms(); ++a2){
	  
	  // no self (nonbonded) interaction
	  if (a2 == it->i || a2 == it->j) continue;

	  DEBUG(13, "\t\thessian " << it->i << " - " << a2);
	  m_nonbonded->calculate_hessian(topo, conf, sim, it->i, a2, h1);

	  DEBUG(13, "\t\thessian " << it->j << " - " << a2);
	  m_nonbonded->calculate_hessian(topo, conf, sim, it->j, a2, h2);

	  for(unsigned int d1=0; d1<3; ++d1){
	    for(unsigned int d2=0; d2<3; ++d2){
	      u1(d1,d2) -= h1(d1,d2);
	      u2(d1,d2) += h2(d1,d2);
	    }
	  }
	} // all hessians

	for(unsigned int d1=0; d1<3; ++d1){
	  for(unsigned int d2=0; d2<3; ++d2){
	    u1(d1,d2) /= topo.mass()(it->i);
	    u2(d1,d2) /= topo.mass()(it->j);
	  }
	}

	math::Vec f1 =   r - math::product(u1, r_nc) * dk * mu / m_parameter[it->type].K;
	math::Vec f2 =  -r - math::product(u2, r_nc) * dk * mu / m_parameter[it->type].K;
	
	// r_nc is normalized!
	f1 += (r_nc - r / r_dist) * dk * mu / (m_parameter[it->type].K * dt2);
	f2 -= (r_nc - r / r_dist) * dk * mu / (m_parameter[it->type].K * dt2);

	assert(force.size() > k && force[k].size() > it->i && force[k].size() > it->j);

	force[k][it->i] = f1;
	force[k][it->j] = f2;

	DEBUG(8, "\t\tr =  " << math::v2s(r));
	DEBUG(8, "\t\tf1 = " << math::v2s(f1));
	DEBUG(8, "\t\tf2 = " << math::v2s(f2));

      }
      else if (a == it->j){
	// already handled...
	continue;
      }
      else{
	DEBUG(12, "\tforces on all other atoms");
	// constraint forces on all other atoms : Eq 27

	DEBUG(13, "\t\thessian " << a << " - " << it->i);
	math::Matrix h1;
	m_nonbonded->calculate_hessian(topo, conf, sim, a, it->i, h1);

	DEBUG(13, "\t\thessian " << a << " - " << it->i);
	math::Matrix h2;
	m_nonbonded->calculate_hessian(topo, conf, sim, a, it->j, h2);
	
	for(unsigned int d1 = 0; d1 < 3; ++d1){
	  for(unsigned int d2 = 0; d2 < 3; ++d2){
	    h1(d1,d2) /= topo.mass()(it->i);
	    h2(d1,d2) /= topo.mass()(it->j);

	    h1(d1,d2) -= h2(d1,d2);
	  }
	}

	assert(force.size() > k && force[k].size() > a);

	// r_nc is normalized
	force[k][a] = math::product(h1, r_nc)  * (-dk) * mu / m_parameter[it->type].K;
	// force[k][a] = 0.0;

	DEBUG(12, "f(" << a << ") = " << math::v2s(force[k][a]));

      }
    }
  }
}

/**
 * shake solute
 */
template<math::boundary_enum B, math::virial_enum V>
void algorithm::Flexible_Constraint::solute
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::vector<std::vector<math::Vec> > & force,
 int & error
 )
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...

  DEBUG(8, "\tFlexible shaking SOLUTE");

  const double start = util::now();

  math::Periodicity<B> periodicity(conf.current().box);
  
  std::vector<double> flex_len;
  calc_distance<B, V>(topo, conf, sim, flex_len);
  // exact solution

  if (sim.param().constraint.solute.flexshake_mode == 2 ||
      sim.param().constraint.solute.flexshake_mode == 3){
    
    calc_undetermined_forces<B, V>(topo, conf, sim, 
				   flex_len, force);
  }
  
  // conf.constraint_force() = 0.0;
  // m_lambda = 0.0;

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);
  
  bool convergence = false;
  while(!convergence){
    DEBUG(9, "\titeration" << std::setw(10) << num_iterations);

    //--------------------------------------------------------------------------------
    // do an iteration
    //--------------------------------------------------------------------------------
    if (sim.param().constraint.solute.flexshake_mode == 0 ||
	sim.param().constraint.solute.flexshake_mode == 1){
      // standard algorithm
      if(iteration<B, V>
	 (topo, conf, convergence, flex_len, skip_now, skip_next,
	  sim.time_step_size(),
	  periodicity, true)
	 ){
	io::messages.add("Flexible SHAKE error. vectors orthogonal",
			 "Flexible_Constraint::solute",
			 io::message::error);
	std::cout << "Flexible SHAKE failure in solute!" << std::endl;
	error = E_SHAKE_FAILURE_SOLUTE;
	return;
      }
    }
    else{
      // exact algorithm
      if(exact_iteration<B, V>
	 (topo, conf, convergence, flex_len, force,
	  sim.time_step_size(), periodicity)
	 ){
	io::messages.add("Exact Flexible SHAKE error. vectors orthogonal",
			 "Flexible_Constraint::solute",
			 io::message::error);
	std::cout << "Exact Flexible SHAKE failure in solute!" << std::endl;
	error =  E_SHAKE_FAILURE_SOLUTE;
	return;
      }
    }
    
    if(++num_iterations > m_max_iterations){
      io::messages.add("Flexible SHAKE error. too many iterations",
		       "Flexible_Constraint::solute",
		       io::message::error);
      error =  E_SHAKE_FAILURE_SOLUTE;
      return;
    }

    skip_now = skip_next;
    skip_next.assign(skip_next.size(), true);

  } // convergence?

  // constraint_force

  // shaken velocity
  for(unsigned int i=0; i<topo.num_atoms(); ++i)
    conf.current().vel(i) = (conf.current().pos(i) - conf.old().pos(i)) / 
      sim.time_step_size();
  
  // calculate the approximate work done
  approx_work<B,V>(topo, conf, flex_len, sim.time_step_size());

  m_timing += util::now() - start;

  error = 0;
}


int algorithm::Flexible_Constraint
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       bool quiet)
{
  if (!quiet){
    std::cout << "FLEXIBLESHAKE\n"
	      << "\tsolute\t";
    if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
      std::cout << "ON\n";
      std::cout << "\t\ttolerance = " 
		<< sim.param().constraint.solute.shake_tolerance << "\n";
      if (sim.param().constraint.solute.flexshake_readin)
	std::cout << "\t\treading velocities along constraints from file\n";
      if (sim.param().constraint.solute.flexshake_mode == 2 ||
	  sim.param().constraint.solute.flexshake_mode == 3)
	std::cout << "\t\tusing the exact algorithm\n";
      else
	std::cout << "\t\tusing the approximate algorithm\n";
      if (sim.param().constraint.solute.flexshake_mode == 0 ||
	  sim.param().constraint.solute.flexshake_mode == 2)
	std::cout << "\t\tusing potential and kinetic energy\n";
      else
	std::cout << "\t\tusing potentialenergy only\n";
      
    }
    else std::cout << "OFF\n";
  
    std::cout << "\tsolvent\t";
  }
  
  if (sim.param().constraint.solvent.algorithm == simulation::constr_flexshake){
    if (!quiet)
      std::cout << "not supported!\n";
    io::messages.add("flexible shake for solvent not implemented", "Flexible_Constraint",
		     io::message::error);
  }
  else if (!quiet) std::cout << "OFF\n";

  if (!quiet)
    std::cout << "END\n";
  
  if (sim.param().start.shake_pos){
    if (!quiet)
      std::cout << "(flexible) shaking initial positions\n";

    // old and current pos and vel are the same...
    conf.old().pos = conf.current().pos;
    conf.old().vel = conf.current().vel;

    // shake the current ones
    if (apply(topo, conf, sim))
      return E_SHAKE_FAILURE;

    // restore the velocities
    conf.current().vel = conf.old().vel;
    
    conf.special().flexible_constraint.
      flexible_vel.assign(conf.special().flexible_constraint.flexible_vel.size(), 0.0);

    // take a step back
    conf.old().pos = conf.current().pos;
    
    if (sim.param().start.shake_vel){
      if (!quiet)
	std::cout << "shaking initial velocities\n";

      for(unsigned int i=0; i<topo.num_atoms(); ++i)
	conf.current().pos(i) = conf.old().pos(i) - 
	  sim.time_step_size() * conf.old().vel(i);
    
      // shake again
      if (apply(topo, conf, sim))
	return E_SHAKE_FAILURE;
    
      // restore the positions
      conf.current().pos = conf.old().pos;
    
      // velocities are in opposite direction (in time)
      for(unsigned int i=0; i<topo.num_atoms(); ++i)
	conf.current().vel(i) = -1.0 * conf.current().vel(i);
      conf.old().vel = conf.current().vel;

      // also change the velocities along the constraints...
      for(size_t i=0; i<conf.special().flexible_constraint.flexible_vel.size(); ++i)
	conf.special().flexible_constraint.flexible_vel[i] *= -1.0;
    }
    
  }
  else if (sim.param().start.shake_vel){
    io::messages.add("shaking velocities without shaking positions illegal.",
		     "shake", io::message::error);
  }

  std::cout << "END\n";
  
  return 0;
}


/**
 * approximate work wasted
 */      
template<math::boundary_enum B, math::virial_enum V>
int algorithm::Flexible_Constraint::approx_work
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 std::vector<double> const & flex_len,
 double dt
 )
{

  if (!conf.special().flexible_constraint.flex_len.size())
    conf.special().flexible_constraint.flex_len = flex_len;

  double w = 0;
  double aw = 0;
  double anw = 0;
  double anaw = 0;
  
  math::Periodicity<B> periodicity(conf.current().box);

  // and constraints
  int k=0;
  for(typename std::vector<topology::two_body_term_struct>::const_iterator
	it = topo.solute().distance_constraints().begin(),
	to = topo.solute().distance_constraints().end();
      it != to;
      ++it, ++k){
    
    DEBUG(10, "i: " << it->i << " j: " << it->j);

    // the position
    math::Vec &ref_i = conf.old().pos(it->i);
    math::Vec &ref_j = conf.old().pos(it->j);

    math::Vec &pos_i = conf.current().pos(it->i);
    math::Vec &pos_j = conf.current().pos(it->j);

    math::Vec dv = 
      (conf.current().vel(it->i) - conf.old().vel(it->i)) -
      (conf.current().vel(it->j) - conf.old().vel(it->j));

    math::Vec r, ref_r;
    periodicity.nearest_image(ref_i, ref_j, ref_r);
    periodicity.nearest_image(pos_i, pos_j, r);

    const double mu = topo.mass()(it->i) * topo.mass()(it->j) /
      (topo.mass()(it->i) + topo.mass()(it->j));

    w += 
      math::dot(dv, ref_r) * mu / dt *
      (sqrt(math::abs2(r)) / sqrt(math::abs2(ref_r)) - 1);

    aw += 
      fabs(math::dot(dv, ref_r) * mu / dt *
	   (sqrt(math::abs2(r)) / sqrt(math::abs2(ref_r)) - 1));

    anw += 
      math::dot(dv, ref_r) * mu / dt *
      (flex_len[k] - conf.special().flexible_constraint.flex_len[k]) /
      sqrt(math::abs2(ref_r));

    anaw += 
      fabs(math::dot(dv, ref_r) * mu / dt *
      (flex_len[k] - conf.special().flexible_constraint.flex_len[k]) /
      sqrt(math::abs2(ref_r)));

  } // constraints

  std::cout.precision(9);
  std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);

  std::cout << "FLEXIBLECONSTRAINT_WORK" << std::setw(30) << w << "\n";
  std::cout << "FLEXIBLECONSTRAINT_ABSWORK" << std::setw(30) << aw << "\n";
  std::cout << "FLEXIBLECONSTRAINT_ANALYTWORK" << std::setw(30) << anw << "\n";
  std::cout << "FLEXIBLECONSTRAINT_ANALYTABSWORK" << std::setw(30) << anaw << "\n";

  // store new lengths...
  conf.special().flexible_constraint.flex_len = flex_len;

  return 0;

}    
