/**
 * @file leap_frog.tcc
 * contains the template methods
 * for the class leap_frog.
 */

/**
 * Leap frog step.
 */
template<typename t_simulation>
void algorithm::leap_frog<t_simulation>
::step(t_simulation &sim,
       interaction::forcefield<t_simulation> &ff,
       double const dt)
{
  // one force calculation per step
  sim.system().exchange_force();
  ff.calculate_interactions(sim);
  
  velocities(sim.system(), sim.topology(), dt);
  positions(sim.system(), dt);
  
}

/**
 * Leap frog velocities
 */
template<typename t_simulation>
void algorithm::leap_frog<t_simulation>
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
template<typename t_simulation>
void algorithm::leap_frog<t_simulation>
::positions(typename t_simulation::system_type &sys,
	    double const dt)
{
  sys.exchange_pos();
  // r = r + v*dt
  sys.pos() = sys.old_pos() + sys.vel() * dt;
  
}
