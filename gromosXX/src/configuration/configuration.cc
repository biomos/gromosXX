/**
 * @file configuration.cc
 * methods definition
 */

#include <stdheader.h>

#include <configuration/configuration_global.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <configuration/configuration.h>

#include <simulation/multibath.h>
#include <simulation/parameter.h>

#include <math/periodicity.h>

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

double configuration_ver = 0.10;

/**
 * Constructor
 */
configuration::Configuration::Configuration()
{
  m_current = &m_state1;
  m_old = &m_state2;

  for(int i=0; i<3; ++i){
    for(int j=0; j<3; ++j){
      current().virial_tensor(i,j) = 0.0;
      old().virial_tensor(i,j) = 0.0;

      current().kinetic_energy_tensor(i,j) = 0.0;
      old().kinetic_energy_tensor(i,j) = 0.0;

      current().pressure_tensor(i,j) = 0.0;
      old().pressure_tensor(i,j) = 0.0;
    }
  }
}

/**
 * copy constructor
 */
configuration::Configuration::Configuration
(
 configuration::Configuration const & conf
 )
{
  m_current = &m_state1;
  m_old = &m_state2;

  for(int i=0; i<3; ++i){
    for(int j=0; j<3; ++j){
      current().virial_tensor(i,j) =
	conf.current().virial_tensor(i,j);
      old().virial_tensor(i,j) =
	conf.old().virial_tensor(i,j);

      current().kinetic_energy_tensor(i,j) = 
	conf.current().kinetic_energy_tensor(i,j);
      old().kinetic_energy_tensor(i,j) = 
	conf.old().kinetic_energy_tensor(i,j);

      current().pressure_tensor(i,j) = 
	conf.current().pressure_tensor(i,j);
      old().pressure_tensor(i,j) = 
	conf.old().pressure_tensor(i,j);
    }
  }
  
  current().pos = conf.current().pos;
  old().pos = conf.old().pos;
  current().posV = conf.current().posV;
  old().posV = conf.old().posV;
  current().vel = conf.current().vel;
  old().vel = conf.old().vel;
  current().force = conf.current().force;
  old().force = conf.old().force;
  current().stochastic_integral = conf.current().stochastic_integral;
  old().stochastic_integral = conf.old().stochastic_integral;
  current().stochastic_seed = conf.current().stochastic_seed;
  old().stochastic_seed = conf.old().stochastic_seed;
  
  current().box = conf.current().box;
  old().box = conf.old().box;
  
  current().energies = conf.current().energies;
  old().energies = conf.old().energies;
  current().averages = conf.current().averages;
  old().averages = conf.old().averages;
  
  current().perturbed_energy_derivatives =
    conf.current().perturbed_energy_derivatives;
  old().perturbed_energy_derivatives =
    conf.old().perturbed_energy_derivatives;
  
  special().dihedral_angle_minimum = conf.special().dihedral_angle_minimum;
  special().flexible_constraint = conf.special().flexible_constraint;
  
  special().jvalue_av = conf.special().jvalue_av;
  special().jvalue_curr = conf.special().jvalue_curr;
  special().jvalue_epsilon = conf.special().jvalue_epsilon;
  
  special().distanceres_av = conf.special().distanceres_av;

  special().pscale = conf.special().pscale;
  
  special().rottrans_constr = conf.special().rottrans_constr;

  special().ramd = conf.special().ramd;
  
  boundary_type = conf.boundary_type;
}

/**
 * operator equal
 */
configuration::Configuration & configuration::Configuration::operator=
(
 configuration::Configuration const & conf
 )
{
  m_current = &m_state1;
  m_old = &m_state2;

  for(int i=0; i<3; ++i){
    for(int j=0; j<3; ++j){
      current().virial_tensor(i,j) =
	conf.current().virial_tensor(i,j);
      old().virial_tensor(i,j) =
	conf.old().virial_tensor(i,j);

      current().kinetic_energy_tensor(i,j) = 
	conf.current().kinetic_energy_tensor(i,j);
      old().kinetic_energy_tensor(i,j) = 
	conf.old().kinetic_energy_tensor(i,j);

      current().pressure_tensor(i,j) = 
	conf.current().pressure_tensor(i,j);
      old().pressure_tensor(i,j) = 
	conf.old().pressure_tensor(i,j);
    }
  }
  
  current().pos = conf.current().pos;
  old().pos = conf.old().pos;
  current().posV = conf.current().posV;
  old().posV = conf.old().posV;
  current().vel = conf.current().vel;
  old().vel = conf.old().vel;
  current().force = conf.current().force;
  old().force = conf.old().force;
  current().stochastic_integral = conf.current().stochastic_integral;
  old().stochastic_integral = conf.old().stochastic_integral;
  current().stochastic_seed = conf.current().stochastic_seed;
  old().stochastic_seed = conf.old().stochastic_seed;
  
  current().box = conf.current().box;
  old().box = conf.old().box;
  
  current().energies = conf.current().energies;
  old().energies = conf.old().energies;
  current().averages = conf.current().averages;
  old().averages = conf.old().averages;
  
  current().perturbed_energy_derivatives =
    conf.current().perturbed_energy_derivatives;
  old().perturbed_energy_derivatives =
    conf.old().perturbed_energy_derivatives;
  
  special().dihedral_angle_minimum = conf.special().dihedral_angle_minimum;
  special().flexible_constraint = conf.special().flexible_constraint;
  
  special().jvalue_av = conf.special().jvalue_av;
  special().jvalue_curr = conf.special().jvalue_curr;
  special().jvalue_epsilon = conf.special().jvalue_epsilon;
  
  special().distanceres_av = conf.special().distanceres_av;

  special().pscale = conf.special().pscale;
  
  special().rottrans_constr = conf.special().rottrans_constr;

  special().ramd = conf.special().ramd;
  
  boundary_type = conf.boundary_type;

  return *this;
}

void configuration::Configuration::init(topology::Topology const & topo,
					simulation::Parameter & param,
					bool gather)
{
  // resize the energy arrays
  const unsigned int num = unsigned(topo.energy_groups().size());
  const unsigned int numb = unsigned(param.multibath.multibath.size());

  DEBUG(5, "number of energy groups: " << num 
	<< "\nnumber of baths: " << numb);

  current().energies.resize(num, numb);
  old().energies.resize(num, numb);
  
  current().energies.ewarn(param.ewarn.limit);
  old().energies.ewarn(param.ewarn.limit);

  current().perturbed_energy_derivatives.resize(num, numb);
  old().perturbed_energy_derivatives.resize(num, numb);

  current().averages.resize(topo, *this, param);
  old().averages.resize(topo, *this, param);

  // possibly resize the dihedral angle monitoring array
  // initialize or set to such a value that it is recalculated in
  // the first step, initialization would require an additional function
  // which can be done only after reading of the coordinates. Would be a
  // nicer solution, but also requires the parameters...
  if(param.print.monitor_dihedrals){
    special().dihedral_angle_minimum.resize
      (topo.solute().dihedrals().size(), 4*math::Pi);
  }
  
  if (param.constraint.solute.algorithm == simulation::constr_flexshake &&
      special().flexible_constraint.flexible_vel.size() == 0){

    special().flexible_constraint.flexible_vel.resize(topo.solute().distance_constraints().size() +
				  topo.perturbed_solute().distance_constraints().size());

    special().flexible_constraint.flexible_ekin.resize(numb);
  }

  if(param.ramd.fc!=0.0){
    special().ramd.force_direction = math::Vec(0.0,0.0,0.0);
    // initialize the ta_average to the minimum distance.
    special().ramd.ta_average = param.ramd.ta_min * exp(1);
    
  }
  
  
  // resize the arrays
  // to make scripting easier...
  resize(topo.num_atoms());

  // gather the molecules!
  // check box size

  // mc: bugfix: chargegroups should be gathered
  //             problem if submolecules are set to 1 atom
  //             (for whatever esoteric reasons)

  if (gather){
    switch(boundary_type){
      case math::vacuum:
	break;
      case math::rectangular:
	{
	  math::Periodicity<math::rectangular> periodicity(current().box);
	  // periodicity.gather_molecules_into_box(*this, topo);
	  periodicity.gather_chargegroups(*this, topo);
	  
	  break;
	}
      case math::triclinic:
	{
	  // NO CUTOFF CHECK -- IMPLEMENT!!!
	  math::Periodicity<math::triclinic> periodicity(current().box);
	  // periodicity.gather_molecules_into_box(*this, topo);
	  periodicity.gather_chargegroups(*this, topo);
	  
	  break;
	}
      case math::truncoct:
	{
	  math::Periodicity<math::truncoct> periodicity(current().box);
	  // periodicity.gather_molecules_into_box(*this, topo);
	  periodicity.gather_chargegroups(*this, topo);

	  break;
	}
      default:
	std::cout << "wrong periodic boundary conditions!";
	io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
    }
  }

  // check periodicity
  switch(boundary_type){
    case math::vacuum:
      break;
    case math::rectangular:
      {
	if (current().box(0)(0) <= 2*param.pairlist.cutoff_long ||
	    current().box(1)(1) <= 2*param.pairlist.cutoff_long ||
	    current().box(2)(2) <= 2*param.pairlist.cutoff_long){
	  io::messages.add("box is too small: not twice the cutoff!",
			   "configuration",
			   io::message::error);
	}
	
	break;
      }
    case math::triclinic:
      {
	// NO CUTOFF CHECK -- IMPLEMENT!!!
	break;
      }
    case math::truncoct:
      {
	if (0.5 * sqrt(3.0) * current().box(0)(0) <= 2 * param.pairlist.cutoff_long){
	  
	  io::messages.add("box is too small: not 4 / sqrt(3) * cutoff!",
			   "configuration",
			   io::message::error);
	}
	break;
      }
    default:
      std::cout << "wrong periodic boundary conditions!";
      io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
  }

  if (boundary_type != math::vacuum){
    if (param.centreofmass.remove_rot){
      io::messages.add("disabling removing of centre of mass rotation (PBC)",
		       "configuration",
		       io::message::notice);
      param.centreofmass.remove_rot = false;
    }
  }
}



/**
 * set the number of atoms.
 */
void configuration::Configuration::resize(unsigned int s)
{
  DEBUG(7, "Configuration resize: " << s);
  
  current().resize(s);
  old().resize(s);
}

/**
 * set the number of atoms.
 * using resizeAndPreserve. Therefore
 * you can enlarge the system (or shrink it)
 * while keeping all existing positions/velocities/...
 * a faster version would be just resize, but then
 * the arrays contain garbage...
 * the energies have to be sized seperately!
 */
void configuration::Configuration::state_struct::resize(unsigned int s)
{
  DEBUG(7, "state struct resize: " << s);

  pos.resize(s);
  posV.resize(s);
  vel.resize(s);
  force.resize(s);
  stochastic_integral.resize(s);
}

namespace configuration
{
  std::ostream &operator<<(std::ostream &os, Configuration &conf)
  {
    os << "a configuration";
    return os;
  }
}

