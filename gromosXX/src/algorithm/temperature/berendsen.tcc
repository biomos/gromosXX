/**
 * @file temperature/berendsen.tcc
 * methods of the berendsen thermostat
 */

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

  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    it = sim.multibath().begin(),
    to = sim.multibath().end();
  
  size_t last = 0;
  
  for( ; it != to; ++it){
    // get a range
    // math::Range bath_range(last, it->last_atom);

    // temperature coupling?
    if (it->tau != -1){

      // calculate the kinetic energy from the old velocities
      double ekin = 0.0;
      for(size_t i=last; i<=it->last_atom; ++i){
	ekin += mass(i) * math::dot(old_vel(i), old_vel(i));
      }
      
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
  math::VArray &vel = sim.system().vel();
  math::VArray const & old_vel = sim.system().old_vel();
  math::SArray const & mass = sim.topology().mass();

  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    it = sim.multibath().begin(),
    to = sim.multibath().end();
  
  size_t last = 0;
  
  for( ; it != to; ++it){
    it->kinetic_energy = 0.0;
    for(size_t i=last; i<=it->last_atom; ++i){
      math::Vec v = 0.5 * (vel(i) + old_vel(i));
      
      it->kinetic_energy += mass(i) * math::dot(v, v);
    } // atoms in bath
    
    it->kinetic_energy *= 0.5;
    
    last = it->last_atom + 1;
    
  } // baths
    
}
  
