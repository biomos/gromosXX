/**
 * @file temperature/berendsen.tcc
 * methods of the berendsen thermostat
 */

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

#include "../../debug.h"

algorithm::Berendsen_Thermostat::Berendsen_Thermostat()
{
}

template<typename t_simulation>
inline void algorithm::Berendsen_Thermostat
::apply(t_simulation &sim, double const dt)
{
  math::VArray &vel = sim.system().vel();
  math::VArray const & old_vel = sim.system().old_vel();
  math::SArray const & mass = sim.topology().mass();

  // calculate the old kinetic energies
  sim.calculate_mol_ekin(false);

  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    it = sim.multibath().begin(),
    to = sim.multibath().end();
  
  size_t last = 0;
  
  for(size_t num=0; it != to; ++it, ++num){
    // get a range
    // math::Range bath_range(last, it->last_atom);

    // temperature coupling?
    if (it->tau != -1){

      // calculate the kinetic energy from the old velocities
      double ekin = 0.0;
      for(size_t i=last; i<=it->last_atom; ++i){
	ekin += mass(i) * math::dot(old_vel(i), old_vel(i));
      }

      assert(sim.system().energies().kinetic_energy.size() > num);
      
      DEBUG(7, "pre-scale ekin: " << 0.5*ekin 
	    << "\tnew: " << sim.system().energies().kinetic_energy[num]);
      
      double free_temp = ekin / (it->dof * math::k_Boltzmann);
      double scale = sqrt(1.0 + dt / it->tau *
		      (it->temperature / free_temp -1));

      // do not calculate the kinetic energy here
      // because SHAKE will be called later on...
      // it->kinetic_energy = 0.0;
      for(size_t i=last; i<=it->last_atom; ++i){
	vel(i) *= scale;
	// math::Vec v = 0.5 * (vel(i) + old_vel(i));
	// it->kinetic_energy += mass(i) * math::dot(v, v);
      }
      
      // it->kinetic_energy *= 0.5;

    }
    // otherwise just do nothing, kinetic energy is calculated later
    // (after SHAKE)

    last = it->last_atom + 1;

  }
  
}

template<typename t_simulation>
inline void algorithm::Berendsen_Thermostat
::calculate_kinetic_energy(t_simulation &sim)
{
  // get rid of previous (temperature scaling) energies
  sim.system().energies().zero(false, true);
  sim.calculate_mol_ekin();
  
  /*
  math::VArray &vel = sim.system().vel();
  math::VArray const & old_vel = sim.system().old_vel();
  math::SArray const & mass = sim.topology().mass();

  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    it = sim.multibath().begin(),
    to = sim.multibath().end();
  
  size_t last = 0;
  
  for(int bath=0; it != to; ++it, ++bath){
    it->kinetic_energy = 0.0;
    for(size_t i=last; i<=it->last_atom; ++i){
      math::Vec v = 0.5 * (vel(i) + old_vel(i));
      
      it->kinetic_energy += mass(i) * math::dot(v, v);
    } // atoms in bath
    
    it->kinetic_energy *= 0.5;
    sim.system().energies().kinetic_energy[bath] = it->kinetic_energy;
    
    last = it->last_atom + 1;
    
  } // baths
  */
}
  
template<typename t_simulation>
inline void algorithm::Berendsen_Thermostat
::calculate_kinetic_energy_lambda_derivative(t_simulation &sim)
{
  math::VArray &vel = sim.system().vel();
  math::VArray const & old_vel = sim.system().old_vel();
  math::SArray const & mass = sim.topology().mass();

  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    it = sim.multibath().begin(),
    to = sim.multibath().end();
  
  size_t last = 0;
  std::vector<double> &e_kin = sim.system().lambda_energies().kinetic_energy;
  
  e_kin.resize(sim.system().energies().kinetic_energy.size());

  for(int bath=0; it != to; ++it, ++bath){
    // assert(e_kin.size() > bath);
    e_kin[bath] = 0.0;
    for(size_t i=last; i<=it->last_atom; ++i){
      
      if(sim.topology().perturbed_atom()[i]){
	DEBUG(7, "\tpertrubed kinetic energy for " << i << " in bath " << bath);
	DEBUG(7, "\tA_mass: " << sim.topology().perturbed_solute().atoms()[i].A_mass() << " B_mass: " << sim.topology().perturbed_solute().atoms()[i].B_mass());
	// for some reason we take the new velocities here
	e_kin[bath] -=
	  (sim.topology().perturbed_solute().atoms()[i].B_mass() -
	   sim.topology().perturbed_solute().atoms()[i].A_mass()) 
	  * math::dot(vel(i), vel(i));
	DEBUG(7, "\tdE_kin/dl: " << e_kin[bath]);
	
      }
      
    } // atoms in bath
    
    e_kin[bath] *= 0.5;
    last = it->last_atom + 1;
    
  } // baths    
}
