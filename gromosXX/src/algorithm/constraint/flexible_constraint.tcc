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
algorithm::Flexible_Constraint<t_simulation>
::Flexible_Constraint(double const tolerance, int const max_iterations)
  : algorithm::Shake<t_simulation>(tolerance, max_iterations),
    m_lfcon(0),
    m_bath(NULL)
{
}

/**
 * shake solute
 */
template<typename t_simulation>
int algorithm::Flexible_Constraint<t_simulation>
::solute(typename simulation_type::topology_type &topo,
	 typename simulation_type::system_type &sys,
	 double dt)
{
  // calculate the flexible constrained distances
  const int first = 0;
  if (m_lfcon)
    calc_distance(topo, sys, first, topo.solute().distance_constraints(), dt);

  // and shake
  return Shake<t_simulation>::solute(topo, sys, dt);
}

template<typename t_simulation>
inline std::vector<double> const &
algorithm::Flexible_Constraint<t_simulation>::r0()const
{
  return m_r0;
}
template<typename t_simulation>
inline std::vector<double> &
algorithm::Flexible_Constraint<t_simulation>::r0()
{
  return m_r0;
}
template<typename t_simulation>
inline std::vector<double> const &
algorithm::Flexible_Constraint<t_simulation>::vel()const
{
  return m_vel;
}
template<typename t_simulation>
inline std::vector<double> &
algorithm::Flexible_Constraint<t_simulation>::vel()
{
  return m_vel;
}
template<typename t_simulation>
inline std::vector<double> const &
algorithm::Flexible_Constraint<t_simulation>::K()const
{
  return m_K;
}
template<typename t_simulation>
inline std::vector<double> &
algorithm::Flexible_Constraint<t_simulation>::K()
{
  return m_K;
}


template<typename t_simulation>
template<typename t_distance_struct>
void algorithm::Flexible_Constraint<t_simulation>
::calc_distance(typename simulation_type::topology_type const &topo,
		typename simulation_type::system_type &sys,
		int const first,
		std::vector<t_distance_struct>
		& constr, double const dt, size_t offset)
{
     
  unsigned int k=offset; //index of flex_constraint_distance
   
  DEBUG(10, "offset: " << k);
  DEBUG(10, "sizes:\n\tvel: " << m_vel.size()
	<< "\n\tK:   " << m_K.size()
	<< "\n\tF:   " << sys.constraint_force().size());

  size_t com, ir;
  sys.energies().flexible_constraints_ir_kinetic_energy.
    assign(sys.energies().kinetic_energy.size(), 0.0);

  //loop over all constraints
  for(typename std::vector<t_distance_struct>
	::iterator
	it = constr.begin(),
	to = constr.end();
      it != to;
      ++it, ++k){
    
    // the position
    assert(sys.pos().size() > first+it->i);
    assert(sys.pos().size() > first+it->j);
    
    math::Vec &pos_i = sys.pos()(first+it->i);
    math::Vec &pos_j = sys.pos()(first+it->j);
	
    math::Vec r;
    sys.periodicity().nearest_image(pos_i, pos_j, r);
    double dist2 = dot(r, r);// actual bond length at (t+Dt) 
    
    assert(sys.old_pos().size() > first+it->i);
    assert(sys.old_pos().size() > first+it->j);

    const math::Vec &ref_i = sys.old_pos()(first+it->i);
    const math::Vec &ref_j = sys.old_pos()(first+it->j);
      
    math::Vec ref_r;
    sys.periodicity().nearest_image(ref_i, ref_j, ref_r);

    double ref_dist2 = dot(ref_r, ref_r);// reference bond length at t 

    // standard formula with velocity along contsraints correction
    // (not the velocityless formula):
    assert(topo.mass().size() > first+it->i);
    assert(topo.mass().size() > first+it->j);
    double red_mass = 1 / (1/topo.mass()(first+it->i) + 1/topo.mass()(first+it->j));
    double dt2 = dt * dt;
      
    // calculate the force on constraint k
    assert(m_vel.size() > k);
    const double force_on_constraint  = (red_mass / dt2) * 
      (sqrt(dist2) - sqrt(ref_dist2) - m_vel[k] * dt);

    // zero energy distance

    // =================================================
    // it->b0: actual flexible constraint distance
    // m_r0:   ideal bond length (no outside force)
    // =================================================

    // double constr_length2 = it->b0 * it->b0;
    const double constr_length2 = m_r0[k] * m_r0[k];
 
     // calculate the flexible constraint distance
    assert(m_K.size() > k);
    it->b0 = force_on_constraint / m_K[k] + sqrt(constr_length2);

    // update the velocity array
    m_vel[k] = (it->b0 - sqrt(ref_dist2)) / dt;

    // now we have to store the kinetic energy of the constraint
    // length change in the correct temperature bath...
    assert(m_bath);
    m_bath->in_bath(it->i, com, ir);
    sys.energies().flexible_constraints_ir_kinetic_energy[ir] +=
      0.5 * red_mass * m_vel[k] * m_vel[k];

    DEBUG(8, "constr " << it->i << " bath " << ir << " ekin " 
	  << 0.5 * red_mass * m_vel[k] * m_vel[k]
	  << " tot ekin " << sys.energies().flexible_constraints_ir_kinetic_energy[ir]);

    // calculate Ekin along the constraints
    Ekin += 0.5 * red_mass * m_vel[k] * m_vel[k];

    // calculate Epot in the bond length constraints
    Epot += 0.5 * m_K[k] * pow(it->b0 - sqrt(constr_length2),2);
      
    // sys.flex_constraint_distance()(k) *= sys.flex_constraint_distance()(k);
    DEBUG(5, "flex_constraint_distance: " << it->b0) ;
        
    } // loop over all constraints
    
}

/**
 * add all bonds to the solute constraint vector and
 * remove them from the bond vector.
 */
template<typename t_simulation>
inline void 
algorithm::Flexible_Constraint<t_simulation>
::add_bond_length_constraints(typename t_simulation::topology_type &topo)
{
  simulation::Solute & solute = topo.solute();

  std::vector<simulation::Bond> bonds;
  std::vector<simulation::Bond>::iterator it = solute.bonds().begin(),
    to = solute.bonds().end();
  for( ; it != to; ++it){
    solute.add_distance_constraint(it->i, it->j,
				   m_bond_parameter[it->type].r0);
    // store in flexible constraints
    m_K.push_back(m_bond_parameter[it->type].K);
    m_r0.push_back(m_bond_parameter[it->type].r0);
  }
  solute.bonds() = bonds;
}

/**
 * add bonds connecting an atom of type iac to the
 * constraint vector and remove from the bond vector.
 */
template<typename t_simulation>
inline void
algorithm::Flexible_Constraint<t_simulation>
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
      solute.add_distance_constraint(it->i, it->j,
				     m_bond_parameter[it->type].r0);
      // store in flexible constraints
      m_K.push_back(m_bond_parameter[it->type].K);
      m_r0.push_back(m_bond_parameter[it->type].r0);
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
algorithm::Flexible_Constraint<t_simulation>
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
      solute.add_distance_constraint(it->i, it->j,
				     m_bond_parameter[it->type].r0);
      // store in flexible constraints
      m_K.push_back(m_bond_parameter[it->type].K);
      m_r0.push_back(m_bond_parameter[it->type].r0);
    }
    else
      bonds.push_back(simulation::Bond(it->i, it->j, it->type));
  }
  solute.bonds() = bonds;
}

template<typename t_simulation>
inline void
algorithm::Flexible_Constraint<t_simulation>
::init(t_simulation &sim, io::Argument &args, io::InTopology &topo,
       io::InInput &input)
{
  // initialize SHAKE / ??
  DEBUG(7, "md init shake");
  int ntc;
  double tolerance;
  input.read_SHAKE(ntc, tolerance);
  this->tolerance(tolerance);

  input.read_FLEXCON(m_lfcon);
  
  switch(ntc){
    case 1:
      break;
    case 2: 
      if (m_lfcon)
	std::cout << "flexible SHAKE bonds containing hydrogen atoms" << std::endl;
      else std::cout << "SHAKE bonds containing hydrogen atoms" << std::endl;

      // read in the parameter
      topo >> *this;
      add_bond_length_constraints(1.0,
				  sim.topology().mass(),
				  sim.topology());
      break;
    case 3: 
      if (m_lfcon) std::cout << "flexible SHAKE all bonds" << std::endl;
      else std::cout << "SHAKE all bonds" << std::endl;
      // read in parameter
      topo >> *this;
      add_bond_length_constraints(sim.topology());
      break;
    default:
      std::cout << "wrong ntc parameter" << std::endl;
  }

  // is there a FLEXCON file?
  if (m_lfcon && args.count("flexcon") == 1){
    std::ifstream flex_file(args["flexcon"].c_str());
    io::InFlexibleConstraints flexin(flex_file);
    flexin.read_FLEXCON(m_vel, sim.topology());
  }
  else{
    // initialize with zero length and zero velocities
    if(m_lfcon){
    io::messages.add("No flexible constraints initial values"
		     " and velocities specified!", "Flexible_Constraints",
		     io::message::warning);
    }
    
    m_vel.assign(m_K.size(), 0.0);
    // the lengths are initialized...
  }

  // give the constraint force the correct size...
  sim.system().constraint_force().resize(sim.topology().solute().
					 distance_constraints().size());

  m_lambda.resize(sim.topology().solute().
		  distance_constraints().size());


  m_bath = &sim.multibath();

}
