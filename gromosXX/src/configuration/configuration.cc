/**
 * @file configuration.cc
 * methods definition
 */

#include <util/stdheader.h>

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

  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j){
      current().virial_tensor(i,j) = 0.0;
      old().virial_tensor(i,j) = 0.0;

      current().kinetic_energy_tensor(i,j) = 0.0;
      old().kinetic_energy_tensor(i,j) = 0.0;

      current().pressure_tensor(i,j) = 0.0;
      old().pressure_tensor(i,j) = 0.0;
    }
  
}

void configuration::Configuration::initialise(topology::Topology & topo,
					      simulation::Parameter const & param)
{
  // resize the energy arrays
  const size_t num = topo.energy_groups().size();
  const size_t numb = param.multibath.multibath.size();

  DEBUG(5, "number of energy groups: " << num 
	<< "\nnumber of baths: " << numb);

  current().energies.resize(num, numb);
  old().energies.resize(num, numb);
  
  current().perturbed_energy_derivatives.resize(num, numb);
  old().perturbed_energy_derivatives.resize(num, numb);

  current().averages.resize(topo, *this, param);
  old().averages.resize(topo, *this, param);

  // resize some special data
  special().rel_mol_com_pos.resize(topo.num_atoms());

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
      special().flexible_vel.size() == 0){

    special().flexible_vel.resize(topo.solute().distance_constraints().size() +
				  topo.perturbed_solute().distance_constraints().size());

    special().flexible_ekin.resize(numb);
  }
  
  // gather the molecules!
  switch(boundary_type){
    case math::vacuum:
      break;
    case math::rectangular:
      {
	math::Periodicity<math::rectangular> periodicity(current().box);
	periodicity.gather_molecules_into_box(*this, topo);
	break;
      }
    case math::triclinic:
      {
	math::Periodicity<math::triclinic> periodicity(current().box);
	periodicity.gather_molecules_into_box(*this, topo);
	break;
      }
    default:
      std::cout << "wrong periodic boundary conditions!";
      io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
  }


}


/**
 * set the number of atoms.
 */
void configuration::Configuration::resize(size_t s)
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
void configuration::Configuration::state_struct::resize(size_t s)
{
  DEBUG(7, "state struct resize: " << s);

  pos.resizeAndPreserve(s);
  vel.resizeAndPreserve(s);
  force.resizeAndPreserve(s);
}

namespace configuration
{
  std::ostream &operator<<(std::ostream &os, Configuration &conf)
  {
    os << "a configuration";
    return os;
  }
}

