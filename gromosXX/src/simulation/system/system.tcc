/**
 * @file system.tcc
 * inline methods definition
 */

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE system

#include "../../debug.h"

/**
 * Constructor
 */
template<math::boundary_enum b>
inline simulation::System<b>::System()
  : m_position1(0),
    m_position2(0),
    m_velocity1(0),
    m_velocity2(0),
    m_force1(0),
    m_force2(0),
    m_periodicity()
{
  m_pos = &m_position1;
  m_old_pos = &m_position2;
  m_vel = &m_velocity1;
  m_old_vel = &m_velocity2;
  m_force = &m_force1;
  m_old_force = &m_force2;

  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j){
      m_virial(i,j) = 0.0;
      m_molecular_kinetic_energy(i,j) = 0.0;
      m_pressure(i,j) = 0.0;
    }
  
}

/**
 * set the number of atoms.
 * using resizeAndPreserve. Therefore
 * you can enlarge the system (or shrink it)
 * while keeping all existing positions/velocities/...
 * a faster version would be just resize, but then
 * the arrays contain garbage...
 */
template<math::boundary_enum b>
inline void simulation::System<b>::resize(size_t s)
{
  DEBUG(7, "system resize: " << s);
  
  m_position1.resizeAndPreserve(s);
  m_position2.resizeAndPreserve(s);
  m_velocity1.resizeAndPreserve(s);
  m_velocity2.resizeAndPreserve(s);
  m_force1.resizeAndPreserve(s);
  m_force2.resizeAndPreserve(s);
  m_rel_mol_com_pos.resizeAndPreserve(s);
  
  m_box_index.resize(s);
  
}

/**
 * position accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::pos()
{
  return *m_pos;
}
/**
 * const position accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::pos()const
{
  return *m_pos;
}
/**
 * old position accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::old_pos()const
{
  return *m_old_pos;
}

/**
 * velocity accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::vel()
{
  return *m_vel;
}
/**
 * const velocity accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::vel()const
{
  return *m_vel;
}

/**
 * const old velocity accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::old_vel()const
{
  return *m_old_vel;
}

/**
 * old velocity accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::old_vel()
{
  return *m_old_vel;
}

/**
 * force accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::force()const
{
  return *m_force;
}

/**
 * force accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::force()
{
  return *m_force;
}

/**
 * old force accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::old_force()const
{
  return *m_old_force;
}

/**
 * constraint force accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::constraint_force()const
{
  return m_constraint_force;
}

/**
 * constraint force accessor
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>::constraint_force()
{
  return m_constraint_force;
}

/**
 * flex constraint distance(const) accessor
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>::flex_constraint_distance()const
{
  return m_flex_constraint_distance;
}

/**
 * flex constraint distance accessor
 */
template<math::boundary_enum b>
inline math::VArray  & simulation::System<b>::flex_constraint_distance()
{
  return m_flex_constraint_distance;
}



/**
 * exchange positions
 */
template<math::boundary_enum b>
inline void simulation::System<b>::exchange_pos()
{
  math::VArray *dummy = m_pos;
  m_pos = m_old_pos;
  m_old_pos = dummy;
}
/**
 * exchange velocities
 */
template<math::boundary_enum b>
inline void simulation::System<b>::exchange_vel()
{
  math::VArray *dummy = m_vel;
  m_vel = m_old_vel;
  m_old_vel = dummy;
}
/**
 * exchange forces
 */
template<math::boundary_enum b>
inline void simulation::System<b>::exchange_force()
{
  math::VArray *dummy = m_force;
  m_force = m_old_force;
  m_old_force = dummy;
}

/**
 * periodicity accessor.
 */
template<math::boundary_enum b>
inline math::Periodicity<b> & simulation::System<b>
::periodicity()
{
  return m_periodicity;
}

/**
 * const periodicity accessor.
 */
template<math::boundary_enum b>
inline math::Periodicity<b> const & simulation::System<b>
::periodicity()const
{
  return m_periodicity;
}
/**
 * const rel_mol_com_pos accessor.
 */
template<math::boundary_enum b>
inline math::VArray const & simulation::System<b>
::rel_mol_com_pos()const
{
  return m_rel_mol_com_pos;
}

/**
 * rel_mol_com_pos accessor.
 */
template<math::boundary_enum b>
inline math::VArray & simulation::System<b>
::rel_mol_com_pos()
{
  return m_rel_mol_com_pos;
}

/**
 * const virial accessor.
 */
template<math::boundary_enum b>
inline math::Matrix const & simulation::System<b>
::virial()const
{
  return m_virial;
}

/**
 * virial accessor.
 */
template<math::boundary_enum b>
inline math::Matrix & simulation::System<b>
::virial()
{
  return m_virial;
}

/**
 * const molecular kinetic energy accessor.
 */
template<math::boundary_enum b>
inline math::Matrix const & simulation::System<b>
::molecular_kinetic_energy()const
{
  return m_molecular_kinetic_energy;
}

/**
 * molecular kinetic energy accessor.
 */
template<math::boundary_enum b>
inline math::Matrix & simulation::System<b>
::molecular_kinetic_energy()
{
  return m_molecular_kinetic_energy;
}

/**
 * const pressure accessor.
 */
template<math::boundary_enum b>
inline math::Matrix const & simulation::System<b>
::pressure()const
{
  return m_pressure;
}

/**
 * pressure accessor.
 */
template<math::boundary_enum b>
inline math::Matrix & simulation::System<b>
::pressure()
{
  return m_pressure;
}

/**
 * box index accessor.
 */
template<math::boundary_enum b>
inline std::vector<typename simulation::System<b>::index_struct> & 
simulation::System<b>::box_indices()
{
  return m_box_index;
}

template<math::boundary_enum b>
inline typename simulation::System<b>::index_struct &
simulation::System<b>::box_index(size_t i)
{
  assert(i < m_box_index.size());
  return m_box_index[i];
}

/**
 * const energy accessor.
 */
template<math::boundary_enum b>
inline simulation::Energy const & 
simulation::System<b>::energies()const
{
  return m_energy;
}

/**
 * energy accessor.
 */
template<math::boundary_enum b>
inline simulation::Energy & 
simulation::System<b>::energies()
{
  return m_energy;
}

/**
 * const energy averages accessor.
 */
template<math::boundary_enum b>
inline simulation::Energy_Average const & 
simulation::System<b>::energy_averages()const
{
  return m_energy_average;
}

/**
 * energy averages accessor.
 */
template<math::boundary_enum b>
inline simulation::Energy_Average & 
simulation::System<b>::energy_averages()
{
  return m_energy_average;
}

/**
 * lambda derivative of the energy accessor.
 */
template<math::boundary_enum b>
inline simulation::Energy & 
simulation::System<b>::lambda_energies()
{
  return m_lambda_energy;
}

/**
 * const lambda derivative averages accessor.
 */
template<math::boundary_enum b>
inline simulation::Energy_Average const & 
simulation::System<b>::lambda_derivative_averages()const
{
  return m_lambda_derivative_average;
}

/**
 * lambda derivative averages accessor.
 */
template<math::boundary_enum b>
inline simulation::Energy_Average & 
simulation::System<b>::lambda_derivative_averages()
{
  return m_lambda_derivative_average;
}

/**
 * check state
 */
template<math::boundary_enum b>
inline int
simulation::System<b>::check_state()const
{
  int result = 0;
  return result;
}


/**
 * calculate the center of mass and the
 * translational kinetic energy of a group of
 * atoms.
 * @param start begin of a group of atoms.
 * @param end of a group of atoms.
 * @mass the masses of (all) atoms.
 * @com_pos returns the center of mass.
 * @com_e_kin returns the tranlational kinetic energy tensor.
 * 
 * @TODO the gathering of the molecule is hardcoded in here.
 * Maybe this should be changed to a generic implementation.
 * Gathering is done in respect to the previous atom. An idea would
 * be to gather as default with respect to the previous atom but
 * letting the user override this (GATHER block).
 * This does not yield the same answer as Phils approach for all cases
 * but maybe for the practical ones???
 */
template<math::boundary_enum b>
void simulation::System<b>::
center_of_mass(Atom_Iterator start, Atom_Iterator end,
	       math::SArray const &mass, 
	       math::Vec &com_pos, math::Matrix &com_e_kin)
{

  com_pos = 0.0;
  double m;
  double tot_mass = 0.0;

  math::Vec p;
  math::Vec prev;
  math::Vec v = 0.0;

  prev = pos()(*start);

  for( ; start != end; ++start){
	
    assert(unsigned(mass.size()) > *start && 
	   unsigned(pos().size()) > *start);

    m = mass(*start);
    tot_mass += m;
    periodicity().nearest_image(pos()(*start), prev, p);
    com_pos += m * (p + prev);
    v += m * vel()(*start);
    prev += p;
  }
  com_pos /= tot_mass;
      
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      com_e_kin(i,j) = 0.5 * v(i) * v(j) / tot_mass;
  
}

template<math::boundary_enum b>
void simulation::System<b>::
molecular_translational_ekin(Atom_Iterator start, Atom_Iterator end,
			     math::SArray const &mass, 
			     math::Vec &com_v, double &com_e_kin,
			     double &e_kin, int mean)
{

  com_v = 0.0;
  com_e_kin = 0.0;
  e_kin = 0.0;
  
  double m;
  double tot_mass = 0.0;
  math::Vec v;

  DEBUG(9, "mol trans ekin: first atom: " << *start);  

  for( ; start != end; ++start){
	
    assert(unsigned(mass.size()) > *start && 
	   unsigned(vel().size()) > *start &&
	   unsigned(old_vel().size()) > *start);

    m = mass(*start);
    tot_mass += m;

    if (mean == 0)
      v = 0.5 * (vel()(*start) + old_vel()(*start));
    else if (mean == -1)
      v = old_vel()(*start);
    else
      v = vel()(*start);
    
    com_v += m * v;
    e_kin += m * dot(v,v);
  }

  com_v /= tot_mass;
      
  com_e_kin = 0.5 * tot_mass * dot(com_v, com_v);
  e_kin *= 0.5;

  DEBUG(9, "com_e_kin: " << com_e_kin << " e_kin: " << e_kin);
  
}


template<math::boundary_enum b>
void simulation::System<b>::
generate_velocities(double const temp, math::SArray const &mass,
		    unsigned int const seed)
{

  math::NormalUnit<double> n;
  n.seed(seed);

  for(int i=0; i<vel().size(); ++i){
    const double sd = sqrt(k_Boltzmann * temp / mass(i));
    for(int d=0; d<3; ++d){
      old_vel()(i)(d) = sd * n.random();
      vel()(i)(d) = old_vel()(i)(d);
    }
  }
}



namespace simulation
{
  template<math::boundary_enum b>
  inline std::ostream &operator<<(std::ostream &os, System<b> &sys)
  {
    os << "a system";
    return os;
  }
}

