/**
 * @file flexible_constraint.tcc
 * implements flexible constraint algorithm
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
algorithm::Perturbed_Flexible_Constraint<t_simulation>
::Perturbed_Flexible_Constraint(double const tolerance, int const max_iterations)
  : algorithm::Flexible_Constraint<t_simulation>(tolerance, max_iterations)
{
}

/**
 * shake solute
 */
template<typename t_simulation>
int algorithm::Perturbed_Flexible_Constraint<t_simulation>
::solute(typename simulation_type::topology_type &topo,
	 typename simulation_type::system_type &sys,
	 double dt)
{
  // calculate the flexible constrained distances
  const int first = 0;

  if (m_lfcon){
    // flexible constraints
    calc_distance(topo, sys, first, 
		  topo.solute().distance_constraints(), dt);
    
    // with an offset
    calc_distance(topo, sys, first, 
		  topo.perturbed_solute().distance_constraints(), dt,
		  topo.solute().distance_constraints().size());
    
  }
   
  // UNFORTUNATELY WE REPEAT HERE THE WHOLE CODE OF
  // PERTURBED SHAKE. MAYBE ONE SHOULD PULL OUT THE
  // LAMBDA DERIVATIVE CALCULATION AND THEN USE MULTIPLE
  // INHERITANCE (FROM FLEXIBLE CONSTRAINTS AND FROM
  // PERTURBED SHAKE)

  // for now shake the whole solute in one go,
  // not bothering about submolecules...
  
  // check whether we shake
  if ((topo.solute().distance_constraints().size() == 0) &&
      (topo.perturbed_solute().distance_constraints().size() == 0))
    return 0;
  
  sys.constraint_force() = 0.0;
  m_lambda = 0.0;

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);
  
  bool convergence = false;
  bool pert_convergence = false;
  
  while(!(convergence && pert_convergence)){

    // non-perturbed constraints
    DEBUG(7, "non-perturbed-flexible-constraints");
    convergence = _shake(topo, sys, first, skip_now, skip_next,
			 topo.solute().distance_constraints(), true);

    // and perturbed constraints
    DEBUG(7, "perturbed-flexible-constraints");
    pert_convergence = _shake(topo, sys, first, skip_now, skip_next,
			      topo.perturbed_solute().distance_constraints(), 
			      true, 
			      topo.solute().distance_constraints().size());

    if(++num_iterations > max_iterations){
      io::messages.add("SHAKE error. too many iterations",
		       "Shake::solute",
		       io::message::critical);
      throw std::runtime_error("SHAKE failure in solute");
    }

    skip_now = skip_next;
    skip_next.assign(skip_next.size(), true);

  } // convergence?

  m_lambda /= -2*dt*dt;

  // constraint_force and free energy derivatives
  // first the non-perturbed part
  size_t k = 0;
  for( ; k < topo.solute().distance_constraints().size(); ++k){

    assert(sys.constraint_force().size() > int(k));
    
    sys.constraint_force()(k) /=  (dt * dt);
    DEBUG(5, "flexible_constraint_force " << 
	  sqrt(dot(sys.constraint_force()(k),
		   sys.constraint_force()(k)) ));

  }

  // now the perturbed part!
  std::vector<simulation::Perturbed_Solute::
    perturbed_distance_constraint_struct>::const_iterator
    it = topo.perturbed_solute().distance_constraints().begin(),
    to = topo.perturbed_solute().distance_constraints().end();


  double m1, m2, dm1, dm2, mu;
  
  for(; it != to; ++it, ++k){
    
    assert(sys.constraint_force().size() > int(k));
    
    sys.constraint_force()(k) /= (dt * dt);
    DEBUG(5, "perturbed flexible_constraint_force " << 
	  sqrt(dot(sys.constraint_force()(k),
		   sys.constraint_force()(k)) ));

    // and the lambda derivative
    assert(topo.mass().size() > it->i + first);
    assert(topo.mass().size() > it->j + first);
    
    m1 = topo.mass()(it->i+first);
    m2 = topo.mass()(it->j+first);
    mu = (m1*m2)/(m1+m2);
    
    assert(int(topo.perturbed_atom().size()) > it->i+first);

    if(topo.perturbed_atom()[it->i+first]){
      dm1 = topo.perturbed_solute().atoms()[it->i+first].B_mass() -
	topo.perturbed_solute().atoms()[it->i+first].A_mass();
    }
    else dm1 = 0;

    if(topo.perturbed_atom()[it->j+first]){
      dm2 = topo.perturbed_solute().atoms()[it->j+first].B_mass() -
	topo.perturbed_solute().atoms()[it->j+first].A_mass();
    }
    else dm2 = 0;

    sys.lambda_energies().constraint_energy +=
      2 * m_lambda(k) * it->b0  *
      (it->B_b0 - it->A_b0 +( it->b0 - m_r0[k] ) *
       ((dm1 + dm2) / (m1 * m2 * mu) - 
	dm2 / m2 - dm1 / m1 - 
	(m_B_K[k] - m_A_K[k]) / m_K[k]
	)
       );

  }
  
  // shaken velocity
  sys.vel() = (sys.pos() - sys.old_pos()) / dt;

  DEBUG(8, "perturbed shake iterations " << num_iterations);

  return num_iterations;

}

/*
template<typename t_simulation>
inline std::vector<double> const &
algorithm::Perturbed_Flexible_Constraint<t_simulation>::A_r0()const
{
  return m_A_r0;
}
template<typename t_simulation>
inline std::vector<double> &
algorithm::Perturbed_Flexible_Constraint<t_simulation>::A_r0()
{
  return m_A_r0;
}
template<typename t_simulation>
inline std::vector<double> const &
algorithm::Perturbed_Flexible_Constraint<t_simulation>::B_r0()const
{
  return m_B_r0;
}
template<typename t_simulation>
inline std::vector<double> &
algorithm::Perturbed_Flexible_Constraint<t_simulation>::B_r0()
{
  return m_B_r0;
}
*/
template<typename t_simulation>
inline std::vector<double> const &
algorithm::Perturbed_Flexible_Constraint<t_simulation>::A_K()const
{
  return m_A_K;
}
template<typename t_simulation>
inline std::vector<double> &
algorithm::Perturbed_Flexible_Constraint<t_simulation>::A_K()
{
  return m_A_K;
}

template<typename t_simulation>
inline std::vector<double> const &
algorithm::Perturbed_Flexible_Constraint<t_simulation>::B_K()const
{
  return m_B_K;
}
template<typename t_simulation>
inline std::vector<double> &
algorithm::Perturbed_Flexible_Constraint<t_simulation>::B_K()
{
  return m_B_K;
}


/*
template<typename t_simulation>
void algorithm::Perturbed_Flexible_Constraint<t_simulation>
::calc_distance(typename simulation_type::topology_type const &topo,
		typename simulation_type::system_type &sys,
		int const first,
		std::vector<simulation::compound::distance_constraint_struct>
		& constr, double const dt)
{
     
  unsigned int k=0;//index of flex_constraint_distance
   
  
  //loop over all constraints
  for(std::vector<simulation::compound::distance_constraint_struct>
	::iterator
	it = constr.begin(),
	to = constr.end();
      it != to;
      ++it, ++k){
    
    // the position
    math::Vec &pos_i = sys.pos()(first+it->i);
    math::Vec &pos_j = sys.pos()(first+it->j);
	
    math::Vec r;
    sys.periodicity().nearest_image(pos_i, pos_j, r);
    double dist2 = dot(r, r);// actual bond length at (t+Dt) 
    

    const math::Vec &ref_i = sys.old_pos()(first+it->i);
    const math::Vec &ref_j = sys.old_pos()(first+it->j);
      
    math::Vec ref_r;
    sys.periodicity().nearest_image(ref_i, ref_j, ref_r);

    double ref_dist2 = dot(ref_r, ref_r);// reference bond length at t 

    // standard formula with velocity along contsraints correction
    // (not the velocityless formula):
    // !!!! still have to check whether the masses are updated !!!!
    double red_mass = 1 / (1/topo.mass()(first+it->i) + 1/topo.mass()(first+it->j));
    double dt2 = dt * dt;
      
    // calculate the force on constraint k
    const double force_on_constraint  = (red_mass / dt2) * 
      (sqrt(dist2) - sqrt(ref_dist2) - m_vel[k] * dt);

    // zero energy distance

    // =================================================
    // it->b0: actual flexible constraint distance
    // m_r0:   ideal bond length (no outside force)
    // =================================================

    // double constr_length2 = it->b0 * it->b0;
    double l = topo.lambda();
    
    const double constr_length2 = (1 - l) * m_A_r0 + l * m_B_r0;
    constr_length2 *= constr_length2;// why first square it and then take the sqrt??????
    
 
     // calculate the flexible constraint distance 
    //coupling param. lambda dependent   
    it->b0 = force_on_constraint /((1 - l) * m_A_K[k] + l * m_B_K[k]) + sqrt(constr_length2);

    // update the velocity array
    m_vel[k] = (it->b0 - sqrt(ref_dist2)) / dt;
    
    // calculate Ekin along the constraints
    Ekin += 0.5 * red_mass * m_vel[k] * m_vel[k];

    // calculate Epot in the bond length constraints
    Epot += 0.5 * ((1 - l) * m_A_K[k] + l * m_B_K[k])  * pow(it->b0 - sqrt(constr_length2),2);
      
    // sys.flex_constraint_distance()(k) *= sys.flex_constraint_distance()(k);
    DEBUG(5, "flex_constraint_distance: " << it->b0) ;
        
    } // loop over all constraints
    
}
*/

/**
 * add all bonds to the solute constraint vector and
 * remove them from the bond vector.
 */
template<typename t_simulation>
inline void 
algorithm::Perturbed_Flexible_Constraint<t_simulation>
::add_bond_length_constraints(typename t_simulation::topology_type &topo)
{

  Flexible_Constraint<t_simulation>::add_bond_length_constraints(topo);
  
  DEBUG(7, "perturbed_flexible_constraints: add_bond_length_constraints");

  simulation::Perturbed_Solute & solute = topo.perturbed_solute();
  
  std::vector<simulation::Perturbed_Bond> bonds;
  std::vector<simulation::Perturbed_Bond>::iterator 
    it = solute.bonds().begin(),
    to = solute.bonds().end();

  for( ; it != to; ++it){
    solute.add_distance_constraint(it->i, it->j,
				   m_bond_parameter[it->type].r0,
				   m_bond_parameter[it->B_type].r0);
    // store in flexible constraints
    m_A_K.push_back(m_bond_parameter[it->type].K);
    m_B_K.push_back(m_bond_parameter[it->B_type].K);

    // m_A_r0.push_back(m_bond_parameter[it->type].A_r0);
    // m_A_r0.push_back(m_bond_parameter[it->type].A_r0);

  }
  solute.bonds().clear();
}


/**
 * add bonds connecting an atom of type iac to the
 * constraint vector and remove from the bond vector.
 */
template<typename t_simulation>
inline void
algorithm::Perturbed_Flexible_Constraint<t_simulation>
::add_bond_length_constraints(int iac,
			      std::vector<int> const &atom_iac,
			      typename t_simulation::topology_type &topo)
{

  Flexible_Constraint<t_simulation>::
    add_bond_length_constraints(iac, atom_iac, topo);
  
  simulation::Perturbed_Solute & solute = topo.perturbed_solute();
  

  std::vector<simulation::Perturbed_Bond> bonds;
  std::vector<simulation::Perturbed_Bond>::iterator
    it = solute.bonds().begin(),
    to = solute.bonds().end();

  for( ; it != to; ++it){
    if(atom_iac[it->i] == iac || atom_iac[it->j] == iac){

      solute.add_distance_constraint(it->i, it->j,
				     m_bond_parameter[it->type].r0,
				     m_bond_parameter[it->B_type].r0);

      // store in flexible constraints
      m_A_K.push_back(m_bond_parameter[it->type].K);
      m_B_K.push_back(m_bond_parameter[it->B_type].K);

      // m_r0.push_back(m_bond_parameter[it->type].r0);
    }
    else
      bonds.push_back(simulation::Perturbed_Bond(*it));
  }
  solute.bonds() = bonds;
}
    
/**
 * add bonds connecting an atom of mass mass to the
 * constraint vector and remove from the bond vector.
 */
template<typename t_simulation>
inline void
algorithm::Perturbed_Flexible_Constraint<t_simulation>
::add_bond_length_constraints(double mass,
			      math::SArray const &atom_mass,
			      typename t_simulation::topology_type &topo)
{
  Flexible_Constraint<t_simulation>::
    add_bond_length_constraints(mass, atom_mass, topo);
  
  simulation::Perturbed_Solute & solute = topo.perturbed_solute();

  std::vector<simulation::Perturbed_Bond> bonds;
  std::vector<simulation::Perturbed_Bond>::iterator 
    it = solute.bonds().begin(),
    to = solute.bonds().end();

  for( ; it != to; ++it){
    if(atom_mass(it->i) == mass || atom_mass(it->j) == mass){

      solute.add_distance_constraint(it->i, it->j,
				     m_bond_parameter[it->type].r0,
				     m_bond_parameter[it->B_type].r0);
      // store in flexible constraints
      m_A_K.push_back(m_bond_parameter[it->type].K);
      m_B_K.push_back(m_bond_parameter[it->B_type].K);

      // m_r0.push_back(m_bond_parameter[it->type].r0);
    }
    else
      bonds.push_back(simulation::Perturbed_Bond(*it));
  }

  solute.bonds() = bonds;
}


template<typename t_simulation>
inline void
algorithm::Perturbed_Flexible_Constraint<t_simulation>
::init(t_simulation &sim, io::Argument &args, io::InTopology &topo,
       io::InInput &input)
{ 
  // initialize SHAKE
  Flexible_Constraint<t_simulation>::init(sim, args, topo, input);
  
  sim.topology().perturbed_solute().
    set_distance_constraints(sim.topology().lambda());

  std::vector<simulation::Perturbed_Solute
    ::perturbed_distance_constraint_struct>
    ::const_iterator
    it = sim.topology().perturbed_solute().distance_constraints().begin(),
    to = sim.topology().perturbed_solute().distance_constraints().end();

  const double l = sim.topology().lambda();
  
  for( ; it != to; ++it){
    m_r0.push_back(it->b0);
  }
  
  for(size_t i=0; i<m_A_K.size(); ++i)
    m_K.push_back((1.0-l)*m_A_K[i]+l*m_B_K[i]);
  
  // give the constraint force the correct size...
  sim.system().constraint_force().resize(sim.topology().solute().
					 distance_constraints().size() + 
					 sim.topology().perturbed_solute().
					 distance_constraints().size());

  // resize only if not read in...
  if (m_vel.size() != sim.topology().solute().distance_constraints().size() +
      sim.topology().perturbed_solute().
      distance_constraints().size()){
    m_vel.assign(sim.topology().solute().distance_constraints().size() +
		 sim.topology().perturbed_solute().
		 distance_constraints().size(), 0.0);
  }

  m_lambda.resize(sim.topology().solute().
		  distance_constraints().size() + 
		  sim.topology().perturbed_solute().
		  distance_constraints().size());

}
