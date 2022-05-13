/**
 * @file perturbed_flexible_constraint.cc
 * contains the template methods for
 * the class Perturbed_Flexible_Constraint.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"
#include "../../interaction/forcefield/forcefield.h"

#include "../../math/periodicity.h"

#include "../../algorithm/constraints/flexible_constraint.h"
#include "../../algorithm/constraints/perturbed_flexible_constraint.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * Constructor.
 */
algorithm::Perturbed_Flexible_Constraint
::Perturbed_Flexible_Constraint
(
 double tolerance,
 int max_iterations,
 interaction::Forcefield *ff
 )
  : Flexible_Constraint(tolerance, max_iterations, ff)
{
}

/**
 * Destructor.
 */
algorithm::Perturbed_Flexible_Constraint
::~Perturbed_Flexible_Constraint()
{
}

/**
 * apply the flexible SHAKE algorithm
 */
int algorithm::Perturbed_Flexible_Constraint
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "applying perturbed flexible SHAKE");
  int error = 0;
  
  conf.special().flexible_constraint.flexible_ekin.assign
    (conf.special().flexible_constraint.flexible_ekin.size(), 0.0);
  
  // check whether we shake
  if (topo.perturbed_solute().distance_constraints().size() && 
      sim.param().constraint.solute.algorithm == 
      simulation::constr_flexshake &&
      sim.param().constraint.ntc > 1){

    DEBUG(8, "\twe need to flexible shake perturbed SOLUTE");

    calc_distance(topo, conf, sim);

    solute(topo, conf, sim, error);
    
    _store_lengths(conf);
    
  }
  
  if (error){
    std::cout << "Perturbed Flexible Constraints: "
	      << "exiting with error condition: E_SHAKE_FAILURE_SOLUTE "
	      << "at step " << sim.steps()
	      << std::endl;
    conf.special().shake_failure_occurred = true;
    return E_SHAKE_FAILURE_SOLUTE;
  }

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

  // return success!
  return 0;
}

/**
 * do one iteration
 */      
template<math::boundary_enum B, math::virial_enum V>
int algorithm::Perturbed_Flexible_Constraint::_iteration
(
 topology::Topology &topo,
 configuration::Configuration & conf,
 bool & convergence,
 std::vector<bool> &skip_now,
 std::vector<bool> &skip_next,
 double dt,
 math::Periodicity<B> const & periodicity
 )
{
  convergence = true;

  // index for constraint_force...
  unsigned int k = 0;
  double const dt2 = dt * dt;
  
  // and constraints
  for(typename std::vector<topology::perturbed_two_body_term_struct>::const_iterator
	it = topo.perturbed_solute().distance_constraints().begin(),
	to = topo.perturbed_solute().distance_constraints().end();
      it != to;
      ++it, ++k ){
	
    // check whether we can skip this constraint
    assert(skip_now.size() > it->i);
    assert(skip_now.size() > it->j);
    
    if (skip_now[it->i] && skip_now[it->j]) continue;

    DEBUG(10, "\ti: " << it->i << " j: " << it->j);

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    // we use the bond lambda
    const double lam = topo.individual_lambda(simulation::bond_lambda)
      [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];
    const double lam_derivative = topo.individual_lambda_derivative
      (simulation::bond_lambda)
      [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];


    // the position
    assert(conf.current().pos.size() > it->i);
    assert(conf.current().pos.size() > it->j);
    
    math::Vec &pos_i = conf.current().pos(it->i);
    math::Vec &pos_j = conf.current().pos(it->j);

    DEBUG(10, "\n\tpositions:");
    DEBUG(10, "\t\ti: " << math::v2s(pos_i) 
	  << "\n\t\tj: " << math::v2s(pos_j));
	
    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);
    double dist2 = math::abs2(r);

    assert(m_perturbed_flex_len.size() > k);

    /*
    if (sqrt(dist2) > 2 * m_perturbed_flex_len[k]){
      std::cerr << "Iteration Error: free distance too large\n"
		<< "\t" << it->i + 1 << " - " << it->j + 1 << "\n"
		<< "\tfree_r = " << sqrt(dist2)
		<< "\tpert_flex_len = " << m_perturbed_flex_len[k]
		<< std::endl;
    }
    */

    DEBUG(10, "\tdist2 = " << dist2);
    DEBUG(10, "\tconstraint " << k << ":");
    DEBUG(10, "\tflex len = " << m_perturbed_flex_len[k]);
    
    double constr_length2 = m_perturbed_flex_len[k] * m_perturbed_flex_len[k];
    double diff = constr_length2 - dist2;

    DEBUG(15, "constr: " << constr_length2 << " dist2: " << dist2);
	  
    if(fabs(diff) >= constr_length2 * m_tolerance * 2.0){
      // we have to shake
      DEBUG(10, "perturbed flexible shaking");
      
      // the reference position
      assert(conf.old().pos.size() > it->i);
      assert(conf.old().pos.size() > it->j);
      
      const math::Vec &ref_i = conf.old().pos(it->i);
      const math::Vec &ref_j = conf.old().pos(it->j);
      
      math::Vec ref_r;
      periodicity.nearest_image(ref_i, ref_j, ref_r);

      /*
      if (math::abs(ref_r) > 2 * sqrt(constr_length2)){
	std::cerr << "Iteration Error: reference distance too large\n"
		  << "\tref_r = " << math::abs(ref_r)
		  << std::endl;
      }
      */

      double sp = dot(ref_r, r);
	  
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
		  << "\tconstr2   : " << constr_length2 << " (" << sqrt(constr_length2) << ")\n"
		  << "\tdist2     : " << dist2 << " (" << sqrt(dist2) << ")\n"
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
      
      // and the lambda derivative
      assert(topo.mass().size() > (it->i));
      assert(topo.mass().size() > (it->j));
    
      const double m1 = topo.mass()(it->i);
      const double m2 = topo.mass()(it->j);
      const double mu = (m1*m2)/(m1+m2);
      double dm1 = 0.0, dm2 = 0.0;

      if(topo.is_perturbed(it->i)){
	dm1 = topo.perturbed_solute().atom(it->i).B_mass() -
	  topo.perturbed_solute().atom(it->i).A_mass();
      }
      else dm1 = 0;

      if(topo.is_perturbed(it->j)){
	dm2 = topo.perturbed_solute().atom(it->j).B_mass() -
	  topo.perturbed_solute().atom(it->j).A_mass();
      }
      else dm2 = 0;

      DEBUG(10, "old dE/dL: " << conf.old().perturbed_energy_derivatives.
	    constraints_energy[topo.atom_energy_group()[it->i]]);
      DEBUG(10, "m1=" << m1 << " m2=" << m2 << " mu=" << mu << " dm1=" 
	    << dm1 << " dm2=" << dm2);

      const double mr0 = (1.0 - lam) * topo.bond_types_harm()[it->A_type].r0 + 
	lam * topo.bond_types_harm()[it->B_type].r0;
      const double mK = (1.0 - lam) * topo.bond_types_harm()[it->A_type].K +
	lam * topo.bond_types_harm()[it->B_type].K;

      DEBUG(10, "mixed r0=" << mr0 << " mixed K=" << mK);

      conf.old().perturbed_energy_derivatives.constraints_energy
	[topo.atom_energy_group()[it->i]] +=
	lam_derivative * (lambda / dt2  * m_perturbed_flex_len[k] *
			  (topo.bond_types_harm()[it->B_type].r0 - 
			   topo.bond_types_harm()[it->A_type].r0 +
			   (m_perturbed_flex_len[k] -  mr0) * 
			   ((dm1 + dm2) / (m1 * m2 * mu) -
			    dm2 / m2 - dm1 / m1 - 
			    (topo.bond_types_harm()[it->B_type].K - 
			     topo.bond_types_harm()[it->A_type].K) / mK
			    )));
	 
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

void algorithm::Perturbed_Flexible_Constraint::calc_distance
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim
 )
{

  Flexible_Constraint::calc_distance(topo, conf, sim);

  SPLIT_BOUNDARY(_calc_distance, topo, conf, sim);
}

template<math::boundary_enum B>
void algorithm::Perturbed_Flexible_Constraint::_calc_distance
(
 topology::Topology const &topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim
 )
{
  DEBUG(8, "\tcalculate perturbed flexible distance");
  
  m_perturbed_flex_len.clear();

  math::Periodicity<B> periodicity(conf.current().box);

  //loop over all constraints
  unsigned int k = unsigned(topo.solute().distance_constraints().size());
  unsigned int com = 0, ir = 0;

  const double dt2 = sim.time_step_size() * sim.time_step_size();
  
  for(std::vector<topology::perturbed_two_body_term_struct>::const_iterator
	it = topo.perturbed_solute().distance_constraints().begin(),
	to = topo.perturbed_solute().distance_constraints().end();
      it != to;
      ++it, ++k){
    
    DEBUG(10, "constraint k=" << k);
    
    // the position
    assert(conf.current().pos.size() > (it->i));
    assert(conf.current().pos.size() > (it->j));
    

    math::Vec const & pos_i = conf.current().pos(it->i);
    math::Vec const & pos_j = conf.current().pos(it->j);
	
    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);

    // unconstrained distance
    const double dist2 = math::abs2(r);
    const double dist = sqrt(dist2);
    
    assert(conf.old().pos.size() > (it->i));
    assert(conf.old().pos.size() > (it->j));

    const math::Vec &ref_i = conf.old().pos(it->i);
    const math::Vec &ref_j = conf.old().pos(it->j);
      
    math::Vec ref_r;
    periodicity.nearest_image(ref_i, ref_j, ref_r);

    // reference distance
    const double ref_dist2 = math::abs2(ref_r);
    const double ref_dist  = sqrt(ref_dist2);
    
    // standard formula with velocity along contsraints correction
    // (not the velocityless formula):
    assert(topo.mass().size() > (it->i));
    assert(topo.mass().size() > (it->j));

    const double red_mass = 1.0 / (1.0/topo.mass()(it->i) + 1.0/topo.mass()(it->j));
      
    // calculate the force on constraint k
    assert(conf.special().flexible_constraint.flexible_vel.size() > k);

    const double force_on_constraint  = (red_mass / dt2) * 
      (dist -  ref_dist - 
       conf.special().flexible_constraint.flexible_vel[k] * sim.time_step_size());
    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    // we use the bond lambda
    const double lam = topo.individual_lambda(simulation::bond_lambda)
      [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];

    /* const double lam_derivative = topo.individual_lambda_derivative
      (simulation::bond_lambda)
      [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]]; */

    /*
    if (fabs(force_on_constraint) > 1E6){
      std::cerr << "ForceOnConstraint Error"
		<< "\n\tf      = " << force_on_constraint
		<< "\n\tdist   = " << dist
		<< "\n\tref_dist = " << ref_dist
		<< "\n\tv        = " << conf.special().flexible_constraint.flexible_vel[k]
		<< "\n\tred_mass = " << red_mass
		<< std::endl;
    }
    */

    // zero energy distance

    // =================================================
    // it->type:  -> ideal bond length
    // flex_len:  flexible constraint distance
    // =================================================

    // const double constr_length2 = topo.bond_types_harm()(it->type).r0 * topo.bond_types_harm()(it->type).r0;
    
    // calculate the flexible constraint distance
    assert(topo.bond_types_harm().size() > it->A_type);
    assert(topo.bond_types_harm().size() > it->B_type);
    
    const double K = (1.0 - lam) * topo.bond_types_harm()[it->A_type].K +
      lam * topo.bond_types_harm()[it->B_type].K;
    const double r0 = (1.0 - lam) * topo.bond_types_harm()[it->A_type].r0 + 
      lam * topo.bond_types_harm()[it->B_type].r0;

    const double new_len = force_on_constraint / K + r0;

    /*
    if (dist > 2 * r0){
      std::cerr << "Calc Distance Error: free dist abnormally large: "
		<< dist
		<< "\tr0 = " << r0
		<< std::endl;
    }

    if (ref_dist > 2 * r0){
      std::cerr << "Calc Distance Error: ref dist abnormally large: "
		<< ref_dist
		<< "\tr0 = " << r0
		<< std::endl;
    }
    
    if (new_len < 0 || new_len > 2 * r0){
      std::cerr << "CALC_DIST ERROR:"
		<< "\n\tnew_len = " << new_len
		<< "\n\tK =       " << K
		<< "\n\tr0 =      " << r0
		<< "\n\tf =       " << force_on_constraint
		<< std::endl;
    }
    */

    // store for shake
    m_perturbed_flex_len.push_back(new_len);

    // update the velocity array
    conf.special().flexible_constraint.flexible_vel[k] = 
      (new_len - ref_dist) / sim.time_step_size();

    // now we have to store the kinetic energy of the constraint
    // length change in the correct temperature bath...
    sim.multibath().in_bath(it->i, com, ir);
    
    assert(conf.special().flexible_constraint.flexible_ekin.size() > ir);
    conf.special().flexible_constraint.flexible_ekin[ir] +=
      0.5 * red_mass * conf.special().flexible_constraint.flexible_vel[k] *
      conf.special().flexible_constraint.flexible_vel[k];
    
    DEBUG(8, "constr " << it->i << " bath " << ir << " ekin " 
	  << 0.5 * red_mass * conf.special().flexible_constraint.flexible_vel[k] * 
	  conf.special().flexible_constraint.flexible_vel[k]);
    
    // calculate Epot in the bond length constraints
    assert(topo.atom_energy_group().size() > it->i);
    assert(conf.old().energies.constraints_energy.size() > topo.atom_energy_group()[it->i]);
    conf.old().energies.constraints_energy[topo.atom_energy_group()[it->i]] += 
      0.5 * K * (r0 - new_len) * (r0 - new_len);
      
    DEBUG(5, "flex_constraint_distance: " << new_len);
    
  } // loop over all constraints
  
}


/**
 * flexible shake perturbed solute
 */
void algorithm::Perturbed_Flexible_Constraint::solute
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 int & error
 )
{
  SPLIT_VIRIAL_BOUNDARY(_solute,
			topo, conf, sim, error);
}

template<math::boundary_enum B, math::virial_enum V>
void algorithm::Perturbed_Flexible_Constraint::_solute
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 int & error
 )
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...

  DEBUG(8, "\tflexible shaking perturbed SOLUTE");
  m_timer.start();

  math::Periodicity<B> periodicity(conf.current().box);
  
  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);
  
  bool convergence = false, pert_convergence = false;
  while(!(convergence && pert_convergence)){

    DEBUG(9, "\titeration" << std::setw(10) << num_iterations);

    if(_iteration<B, V>
       (topo, conf, pert_convergence, skip_now, skip_next,
	sim.time_step_size(), periodicity))
      {
	std::cout << "Perturbed Flexible SHAKE failure in solute!" << std::endl;
	std::cout << "after " << num_iterations << " iterations" << std::endl;

	io::messages.add("Perturbed Flexible SHAKE error. vectors orthogonal",
			 "Perturbed_Flexible_Constraint::solute",
			 io::message::error);
	error = E_SHAKE_FAILURE_SOLUTE;
	return;
      }

    // pert_convergence = true;

    if(Flexible_Constraint::_iteration<B, V>
       (topo, conf, convergence, skip_now, skip_next,
	sim.time_step_size(), periodicity))
      {
	std::cout << "Flexible SHAKE failure in solute!" << std::endl;
	std::cout << "after " << num_iterations << " iterations" << std::endl;

	io::messages.add("Flexible SHAKE error. vectors orthogonal",
			 "(Perturbed) Flexible_Constraint::solute",
			 io::message::error);
	error = E_SHAKE_FAILURE_SOLUTE;
	return;
      }
    // convergence = true;
    
    if(++num_iterations > m_max_iterations){
      io::messages.add("Perturbed Flexible SHAKE error. too many iterations",
		       "Perturbed_Flexible_Constraint::solute",
		       io::message::error);
      error = E_SHAKE_FAILURE_SOLUTE;
      return;
    }

    skip_now = skip_next;
    skip_next.assign(skip_next.size(), true);

  } // convergence?

  m_timer.stop();
  error = 0;
} // solute


int algorithm::Perturbed_Flexible_Constraint
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       std::ostream & os,
       bool quiet)
{
  if (!quiet){
    os << "PERTURBED_FLEXIBLESHAKE\n"
	      << "\tsolute\t";
    if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
      os << "ON\n";
      os << "\t\ttolerance = " << sim.param().constraint.solute.shake_tolerance << "\n";
      if (sim.param().constraint.solute.flexshake_readin)
	os << "\t\treading velocities along constraints from file\n";
    }
    else os << "OFF\n";
    
    os << "\tsolvent\t";
  }
  
  if (sim.param().constraint.solvent.algorithm == simulation::constr_flexshake){
    os << "not supported!\n";
    io::messages.add("flexible shake for solvent not implemented", "Flexible_Constraint",
		     io::message::error);
  }
  else if (!quiet) os << "OFF\n";
  
  if (!quiet)
    os << "END\n";

  if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake) {
    // loop over the constraints to find out which atoms are constrained
    {
      std::vector<topology::two_body_term_struct>::const_iterator
      it = topo.solute().distance_constraints().begin(),
              to = topo.solute().distance_constraints().end();
      for (; it != to; ++it) {
        constrained_atoms().insert(it->i);
        constrained_atoms().insert(it->j);
      }
    }
    // loop over the perturbed constraints
    {
      std::vector<topology::perturbed_two_body_term_struct>::const_iterator
      it = topo.perturbed_solute().distance_constraints().begin(),
              to = topo.perturbed_solute().distance_constraints().end();
      for (; it != to; ++it) {
        constrained_atoms().insert(it->i);
        constrained_atoms().insert(it->j);
      }
    }
  }
  
  if (sim.param().start.shake_pos) {
    if (!quiet)
      os << "\n\tflexible shaking initial positions\n";

    // old and current pos and vel are the same for constrained atoms...
    std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.old().pos(*it) = conf.current().pos(*it);
      conf.old().vel(*it) = conf.current().vel(*it);
    }

    // shake the current ones
    if (apply(topo, conf, sim))
      return E_SHAKE_FAILURE;

    it = constrained_atoms().begin();
    for (; it != to; ++it) {
      // restore the velocities
      conf.current().vel(*it) = conf.old().vel(*it);
      // take a step back
      conf.old().pos(*it) = conf.current().pos(*it);
    }

    if (sim.param().start.shake_vel) {
      if (!quiet)
        os << "\tshaking initial velocities\n";

      it = constrained_atoms().begin();
      for (; it != to; ++it) {
        conf.current().pos(*it) = conf.old().pos(*it) -
                sim.time_step_size() * conf.old().vel(*it);
      }

      // shake again
      if (apply(topo, conf, sim))
        return E_SHAKE_FAILURE;

      it = constrained_atoms().begin();
      for (; it != to; ++it) {
        // restore the positions
        conf.current().pos(*it) = conf.old().pos(*it);
        // velocities are in opposite direction (in time)
        conf.current().vel(*it) = -1.0 * conf.current().vel(*it);
        conf.old().vel(*it) = conf.current().vel(*it);
      }
    } // if shake vel
  } else if (sim.param().start.shake_vel) {
    io::messages.add("shaking velocities without shaking positions illegal.",
            "shake", io::message::error);
  }
  
  if (!quiet)
    os << "END\n";
  
  return 0;
}

void algorithm::Perturbed_Flexible_Constraint::_store_lengths
(
 configuration::Configuration & conf
 )
{

  if (conf.special().flexible_constraint.flex_len.size() < 
      m_perturbed_flex_len.size() + m_flex_len.size())
    conf.special().flexible_constraint.flex_len.resize(m_flex_len.size() + 
						       m_perturbed_flex_len.size());

  unsigned int k=0;
  for( ; k<m_flex_len.size(); ++k){
    assert(conf.special().flexible_constraint.flex_len.size() > k);
    assert(m_flex_len.size() > k);
    conf.special().flexible_constraint.flex_len[k] = m_flex_len[k];
  }
  

  for(unsigned int kk=0; kk<m_perturbed_flex_len.size(); ++k, ++kk){
    assert(conf.special().flexible_constraint.flex_len.size() > k);
    assert(m_perturbed_flex_len.size() > kk);
    conf.special().flexible_constraint.flex_len[k] = m_perturbed_flex_len[kk];
  }
  
}
