/**
 * @file state_properties.cc
 * functions to calculate state properties.
 */

#include <util/stdheader.h>

#include <configuration/configuration_global.h>

#include <algorithm/algorithm.h>
#include <configuration/configuration.h>
#include <topology/topology.h>

#include <math/periodicity.h>

#include <configuration/state_properties.h>

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

template <math::boundary_enum b>
static void
_center_of_mass(topology::Atom_Iterator start, topology::Atom_Iterator end,
		math::SArray const &mass, math::Vec &com_pos, 
		math::Matrix &com_e_kin, configuration::Configuration const & conf)
{
  math::Periodicity<b> periodicity(conf.current().box);
  
  com_pos = 0.0;
  double m;
  double tot_mass = 0.0;

  math::Vec p;
  math::Vec prev;
  math::Vec v = 0.0;

  math::VArray const & pos = conf.current().pos;
  math::VArray const & vel = conf.current().vel;
  prev = pos(*start);
  

  for( ; start != end; ++start){
	
    assert(unsigned(mass.size()) > *start && 
	   unsigned(pos.size()) > *start);

    m = mass(*start);
    tot_mass += m;
    periodicity.nearest_image(pos(*start), prev, p);
    com_pos += m * (p + prev);
    v += m * vel(*start);
    prev += p;
  }

  com_pos /= tot_mass;
      
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      com_e_kin(i,j) = 0.5 * v(i) * v(j) / tot_mass;
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
void configuration::State_Properties::
center_of_mass(topology::Atom_Iterator start, topology::Atom_Iterator end,
	       math::SArray const &mass, 
	       math::Vec &com_pos, math::Matrix &com_e_kin)
{


  switch(m_configuration.boundary_type){
    case math::vacuum :
      _center_of_mass<math::vacuum>(start, end, mass, com_pos, com_e_kin,
				    m_configuration);
      break;
    case math::triclinic :
      _center_of_mass<math::triclinic>(start, end, mass, com_pos, com_e_kin,
				       m_configuration);
      break;
    case math::rectangular :
      _center_of_mass<math::rectangular>(start, end, mass, com_pos, com_e_kin,
					 m_configuration);
      break;
    default:
      throw std::string("Wrong boundary type");
  }
  
}

void configuration::State_Properties::
molecular_translational_ekin(topology::Atom_Iterator start, 
			     topology::Atom_Iterator end,
			     math::SArray const &mass, 
			     math::Vec &com_v, double &com_e_kin,
			     double &e_kin, math::Vec &new_com_v,
			     double &new_com_e_kin, double &new_e_kin)
{

  com_v = 0.0;
  com_e_kin = 0.0;
  e_kin = 0.0;
  new_com_v = 0.0;
  new_com_e_kin = 0.0;
  new_e_kin = 0.0;
  
  double m;
  double tot_mass = 0.0;
  math::Vec v, new_v;

  DEBUG(9, "mol trans ekin: first atom: " << *start);  

  math::VArray const & vel = m_configuration.current().vel;
  math::VArray const & old_vel = m_configuration.old().vel;

  for( ; start != end; ++start){
	
    assert(unsigned(mass.size()) > *start && 
	   unsigned(vel.size()) > *start &&
	   unsigned(old_vel.size()) > *start);

    m = mass(*start);
    tot_mass += m;

    v = 0.5 * (vel(*start) + old_vel(*start));
    new_v = vel(*start);
    
    com_v += m * v;
    e_kin += m * dot(v,v);
    DEBUG(11, "scaling ekin mass=" << m << " v=" << new_v);
    DEBUG(11, "av v="<<v);
    DEBUG(11, "old_v=" << old_vel(*start));
    
    new_com_v += m * new_v;
    new_e_kin += m * dot(new_v, new_v);

  }

  com_v /= tot_mass;
  new_com_v /= tot_mass;
  
  com_e_kin = 0.5 * tot_mass * dot(com_v, com_v);
  e_kin *= 0.5;
  
  new_com_e_kin = 0.5 * tot_mass * math::dot(new_com_v, new_com_v);
  new_e_kin *= 0.5;

  DEBUG(9, "com_e_kin: " << com_e_kin << " e_kin: " << e_kin);
  DEBUG(9, "new com_e_kin: " << new_com_e_kin << " new e_kin: " << new_e_kin);
  
}
