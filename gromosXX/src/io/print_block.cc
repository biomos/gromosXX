/**
 * @file print_block.cc
 * routines to print out the various blocks.
 */

#include "../math/gmath.h"
#include "../io/message.h"

#include "../simulation/simulation.h"

/**
 * Print the multibath block.
 * @TODO write summary information.
 */
std::ostream & operator<<(std::ostream &os, simulation::Multibath &bath)
{
  os << "\nMULTIBATH\n";
  
  os.precision(2);
  os.setf(std::ios_base::fixed, std::ios_base::floatfield);
  
  os << std::setw(10) << "LAST"
     << std::setw( 8) << "TEMP0"
     << std::setw( 6) << "TAU"
     << std::setw(10) << "DOF"
     << std::setw(10) << "SOLUC"
     << std::setw(10) << "SOLVC"
     << std::setw(10) << "EKIN"
     << std::setw(10) << "TEMP"
     << "\n";
  
  std::vector<simulation::bath_struct>::const_iterator
    it = bath.begin(),
    to = bath.end();
  
  for( ; it != to; ++it){
    os << std::setw(10) << it->last_atom + 1
       << std::setw( 8) << it->temperature
       << std::setw( 6) << it->tau
       << std::setw(10) << it->dof
       << std::setw(10) << it->solute_constr_dof
       << std::setw(10) << it->solvent_constr_dof
       << std::setw(10) << it->kinetic_energy;
    if (it->kinetic_energy == 0){
      os << std::setw(10) << 0;
    }
    else{
      os << std::setw(10) 
	 << 2 * it->kinetic_energy / (math::k_Boltzmann * it->dof);
    }
    
    os << "\n";

  }

  os << "END\n";
  return os;
}
