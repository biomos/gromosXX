/**
 * @file leap_frog.tcc
 * contains the template methods
 * for the class leap_frog.
 */

/**
 * Constructor.
 */
template<typename t_simulation, typename t_thermostat>
algorithm::Leap_Frog<t_simulation, t_thermostat>
::Leap_Frog()
  : m_thermostat()
{
}

/**
 * Leap frog step.
 */
template<typename t_simulation, typename t_thermostat>
void algorithm::Leap_Frog<t_simulation, t_thermostat>
::step(t_simulation &sim,
       interaction::Forcefield<t_simulation> &ff,
       double const dt)
{
  // one force calculation per step
  sim.system().exchange_force();
  ff.calculate_interactions(sim);
  
  velocities(sim.system(), sim.topology(), dt);

  m_thermostat.apply(sim, dt);

  positions(sim.system(), dt);
  
}

/**
 * Leap frog velocities
 */
template<typename t_simulation, typename t_thermostat>
void algorithm::Leap_Frog<t_simulation, t_thermostat>
::velocities(typename t_simulation::system_type &sys,
	     typename t_simulation::topology_type &topo,
	     double const dt)
{
  sys.exchange_vel();
  // v = v + f * dt / m
  sys.vel() = sys.old_vel() + sys.force() * dt / topo.mass();
}

/**
 * Leap frog positions
 */
template<typename t_simulation, typename t_thermostat>
void algorithm::Leap_Frog<t_simulation, t_thermostat>
::positions(typename t_simulation::system_type &sys,
	    double const dt)
{
  sys.exchange_pos();
  // r = r + v*dt
  sys.pos() = sys.old_pos() + sys.vel() * dt;
  
}
