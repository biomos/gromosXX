/**
 * @file leap_frog.tcc
 * contains the template methods
 * for the class leap_frog.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration
#include "../../debug.h"

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
template<typename t_interaction_spec>
void algorithm::Leap_Frog<t_simulation, t_thermostat>
::step(t_simulation &sim,
       interaction::Forcefield<t_simulation, t_interaction_spec> &ff,
       double const dt)
{
  // one force calculation per step
  DEBUG(5, "Leap frog step");
  
  sim.system().exchange_force();
  DEBUG(7, "Leap frog: calculate interactions");
  ff.calculate_interactions(sim);
  
  DEBUG(7, "Leap frog: velocities");
  velocities(sim.system(), sim.topology(), dt);

  DEBUG(7, "Leap frog: thermostat");
  m_thermostat.apply(sim, dt);

  DEBUG(7, "Leap frog: positions");
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
