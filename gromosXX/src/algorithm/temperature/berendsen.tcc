/**
 * @file berendsen.tcc
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

      it->kinetic_energy = 0.0;
      for(size_t i=last; i<=it->last_atom; ++i){
	vel(i) *= scale;
	math::Vec v = 0.5 * (vel(i) + old_vel(i));
	it->kinetic_energy += mass(i) * math::dot(v, v);
      }
      
      it->kinetic_energy *= 0.5;

    }
    else{
      // only calculate kinetic energy
      // it is the average over the velocities before and after
      // the old position

      it->kinetic_energy = 0.0;
      for(size_t i=last; i<=it->last_atom; ++i){
	
	math::Vec v = 0.5 * (vel(i) + old_vel(i));

	it->kinetic_energy += mass(i) * math::dot(v, v);
      }
      it->kinetic_energy *= 0.5;

    }

    last = it->last_atom;

  }
  
}
