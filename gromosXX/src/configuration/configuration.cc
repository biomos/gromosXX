/**
 * @file configuration.cc
 * methods definition
 */

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

#include <util/stdheader.h>
#include <configuration/configuration_global.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>

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

