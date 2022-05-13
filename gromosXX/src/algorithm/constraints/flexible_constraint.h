/**
 * @file flexible_constraint.h
 * the flexible shake algorithm
 */

#ifndef INCLUDED_FLEXIBLE_CONSTRAINT_H
#define INCLUDED_FLEXIBLE_CONSTRAINT_H

namespace interaction
{
  class Nonbonded_Interaction;
}

namespace algorithm
{
  /**
   * @class Flexible_Constraint
   * calculates the flexible constraint distance
   */
  class Flexible_Constraint : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Flexible_Constraint(double const tolerance = 0.000001,
			int const max_iterations = 1000,
			interaction::Forcefield *ff = NULL);

    /**
     * Destructor.
     */
    virtual ~Flexible_Constraint();

    /**
     * initialization.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

    /**
     * apply flexible shake.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    
    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

    /**
     * tolerance.
     */
    double const & tolerance()const {return m_tolerance;}
    /**
     * max iterations.
     */
    int const & max_iterations()const {return m_max_iterations;}

    /**
     * accessor to the constrained atoms
     */
    std::set<unsigned int> & constrained_atoms() {
      return m_constrained_atoms;
    }
     /**
     * accessor to the constrained atoms
     */
    const std::set<unsigned int> & constrained_atoms() const {
      return m_constrained_atoms;
    }

    void calc_distance
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim
     );

    void solute
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     int & error
     );

    void _store_lengths
    (
     configuration::Configuration & conf
     );

  protected:
    /**
     * shake tolerance
     */
    double m_tolerance;
    /**
     * max iterations.
     */
    const int m_max_iterations;
    /**
     * the atoms that are involved in the contraints
     */
    std::set<unsigned int> m_constrained_atoms;
    /**
     * flexible constraint lengths
     */
    std::vector<double> m_flex_len;
    
    /**
     * the nonbonded's for exact algorithm...
     * (only implemented for diatomic molecules...
     *  no bonded terms!)
     */
    interaction::Nonbonded_Interaction * m_nonbonded;

    /**
     * the (undetermined) constrained forces (without the lambda factor)
     */
    std::vector<std::vector<math::Vec> > m_force;

    template<math::boundary_enum B, math::virial_enum V>
    int _iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     bool & convergence,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     double dt,
     math::Periodicity<B> const & periodicity
     );

    template<math::boundary_enum B, math::virial_enum V>
    int _exact_iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     bool & convergence,
     double dt,
     math::Periodicity<B> const & periodicity
     );

    template<math::boundary_enum B>
    void _calc_distance
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim
     );

    template<math::boundary_enum B, math::virial_enum V>
    void _calc_undetermined_forces
    (
     topology::Topology &topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     );

    template<math::boundary_enum B, math::virial_enum V>
    void _solute
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     int & error
     );

    template<math::boundary_enum B, math::virial_enum V>
    int _approx_work
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     double dt
     );

  };
  
} // algorithm


/*********************************************************************
 * Definitions of templated functions
 *********************************************************************/

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/nonbonded/interaction/nonbonded_interaction.h"

#include "../../interaction/forcefield/forcefield.h"

#include "../../math/periodicity.h"

#include "../../algorithm/constraints/flexible_constraint.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"

// save previous values
#ifdef MODULE
#define MODULE_PREV MODULE
#endif

#ifdef SUBMODULE
#define SUBMODULE_PREV SUBMODULE
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * do one iteration
 */      
template<math::boundary_enum B, math::virial_enum V>
int algorithm::Flexible_Constraint::_iteration
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 bool & convergence,
 std::vector<bool> & skip_now,
 std::vector<bool> & skip_next,
 double dt,
 math::Periodicity<B> const & periodicity
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
	
    double constr_length2 = m_flex_len[k] * m_flex_len[k];
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
		  << "\tatom i    : " << it->i + 1 << "\n"
		  << "\tatom j    : " << it->j + 1 << "\n"
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

      if (V == math::atomic_virial) {
        for (int a = 0; a < 3; ++a) {
          for (int aa = 0; aa < 3; ++aa) {
            conf.old().virial_tensor(a, aa) +=
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
int algorithm::Flexible_Constraint::_exact_iteration
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 bool & convergence,
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
    
    double constr_length2 = m_flex_len[k] * m_flex_len[k];
    double diff = constr_length2 - dist2;

    DEBUG(12, "constr: " << constr_length2 << " dist2: " << dist2);
    
    if(fabs(diff) >= constr_length2 * m_tolerance * 2.0){
      // we have to shake
      DEBUG(10, "exact flexible shaking");
      
      // the undetermined constraint forces
      // (as function of the reference position...)
      const math::Vec &f_i = m_force[k][it->i];
      const math::Vec &f_j = m_force[k][it->j];
      
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

      if (V == math::atomic_virial) {
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
	conf.current().pos(a) -= 2 * dt2 / topo.mass()(a) * lambda * m_force[k][a];
      }
	  
      convergence = false;
      
    } // we have to shake
  } // constraints
      
  
  return 0;

}

//================================================================================
// flexible constraint distance
//================================================================================

template<math::boundary_enum B>
void algorithm::Flexible_Constraint::_calc_distance
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim
 )
{

  DEBUG(7, "FLEXIBLE CONSTRAINTS::calc_distance");
  
  m_flex_len.clear();
  
  math::Periodicity<B> periodicity(conf.current().box);

  //loop over all constraints
  unsigned int k = 0;
  unsigned int com = 0, ir = 0;
  
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
    double force_on_constraint = 0.0;

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

    // const double constr_length2 = topo.bond_types_harm()(it->type).r0 * topo.bond_types_harm()(it->type).r0;
    
    // calculate the flexible constraint distance
    DEBUG(10, "F(c) = " << force_on_constraint);
    DEBUG(10, "it->type = " << it->type << " of " << topo.bond_types_harm().size());
    
    assert(topo.bond_types_harm().size() > it->type);
    DEBUG(10, "K = " << topo.bond_types_harm()[it->type].K << "     r0 = " << topo.bond_types_harm()[it->type].r0);

    const double new_len = force_on_constraint / topo.bond_types_harm()[it->type].K + 
      topo.bond_types_harm()[it->type].r0;
    
    // store for shake
    m_flex_len.push_back(new_len);

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
      0.5 * topo.bond_types_harm()[it->type].K * (topo.bond_types_harm()[it->type].r0 - new_len) * 
      (topo.bond_types_harm()[it->type].r0 - new_len);
      
    DEBUG(5, "flex_constraint_distance: " << new_len);
    
  } // loop over all constraints
}

//================================================================================
// undetermined forces
//================================================================================
template<math::boundary_enum B, math::virial_enum V>
void algorithm::Flexible_Constraint::_calc_undetermined_forces
(
 topology::Topology &topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  DEBUG(7, "FLEXIBLE CONSTRAINTS: _calc_undetermined_forces");
  
  // resize to number of constraints
  m_force.resize(topo.solute().distance_constraints().size());

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
    
    m_force[k].assign(topo.num_atoms(), math::Vec(0.0));

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
      assert(m_flex_len.size() > k);
      const double dk = m_flex_len[k];
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

	math::Vec f1 =   r - math::product(u1, r_nc) * dk * mu / topo.bond_types_harm()[it->type].K;
	math::Vec f2 =  -r - math::product(u2, r_nc) * dk * mu / topo.bond_types_harm()[it->type].K;
	
	// r_nc is normalized!
	f1 += (r_nc - r / r_dist) * dk * mu / (topo.bond_types_harm()[it->type].K * dt2);
	f2 -= (r_nc - r / r_dist) * dk * mu / (topo.bond_types_harm()[it->type].K * dt2);

	assert(m_force.size() > k && m_force[k].size() > it->i && m_force[k].size() > it->j);

	m_force[k][it->i] = f1;
	m_force[k][it->j] = f2;

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

	assert(m_force.size() > k && m_force[k].size() > a);

	// r_nc is normalized
	m_force[k][a] = math::product(h1, r_nc)  * (-dk) * mu / topo.bond_types_harm()[it->type].K;
	DEBUG(12, "f(" << a << ") = " << math::v2s(m_force[k][a]));

      }
    }
  }
}

template<math::boundary_enum B, math::virial_enum V>
void algorithm::Flexible_Constraint::_solute
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 int & error
 )
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...

  DEBUG(8, "\tFlexible shaking SOLUTE");

  m_timer.start();

  math::Periodicity<B> periodicity(conf.current().box);
  
  // exact solution
  if (sim.param().constraint.solute.flexshake_mode == 2 ||
      sim.param().constraint.solute.flexshake_mode == 3){
    
    _calc_undetermined_forces<B, V>(topo, conf, sim);
  }
  
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
      if(_iteration<B, V>
	 (topo, conf, convergence, skip_now, skip_next,
	  sim.time_step_size(),
	  periodicity)
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
      if(_exact_iteration<B, V>
	 (topo, conf, convergence, sim.time_step_size(), periodicity)
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

  // shaken velocity:
  // stochastic dynamics, energy minimisation, analysis needs to shake without
  // velocity correction (once; it shakes twice...)
  if (!sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze) {
    std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.current().vel(*it) = (conf.current().pos(*it) - conf.old().pos(*it)) /
              sim.time_step_size();
    }
  }
  
  // calculate the approximate work done
  _approx_work<B,V>(topo, conf, sim.time_step_size());

  // store the constraint lengths
  _store_lengths(conf);

  m_timer.stop();

  error = 0;
}

/**
 * approximate work wasted
 */      
template<math::boundary_enum B, math::virial_enum V>
int algorithm::Flexible_Constraint::_approx_work
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 double dt
 )
{

  if (!conf.special().flexible_constraint.flex_len.size())
    conf.special().flexible_constraint.flex_len = m_flex_len;

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
      (m_flex_len[k] - conf.special().flexible_constraint.flex_len[k]) /
      sqrt(math::abs2(ref_r));

    anaw += 
      fabs(math::dot(dv, ref_r) * mu / dt *
      (m_flex_len[k] - conf.special().flexible_constraint.flex_len[k]) /
      sqrt(math::abs2(ref_r)));

  } // constraints

  std::cout.precision(9);
  std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);

  std::cout << "FLEXIBLECONSTRAINT_WORK" << std::setw(30) << w << "\n";
  std::cout << "FLEXIBLECONSTRAINT_ABSWORK" << std::setw(30) << aw << "\n";
  std::cout << "FLEXIBLECONSTRAINT_ANALYTWORK" << std::setw(30) << anw << "\n";
  std::cout << "FLEXIBLECONSTRAINT_ANALYTABSWORK" << std::setw(30) << anaw << "\n";

  // store new lengths...
  // conf.special().flexible_constraint.flex_len = m_flex_len;

  return 0;

}

// set macros to previous values if needed
#ifdef MODULE_PREV
#undef MODULE
#define MODULE MODULE_PREV
#undef MODULE_PREV
#endif

#ifdef SUBMODULE_PREV
#undef SUBMODULE
#define SUBMODULE SUBMODULE_PREV
#undef SUBMODULE_PREV
#endif


#endif
