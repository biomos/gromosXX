/**
 * @file runge_kutta.tcc
 * contains the template methods
 * for the class runge_kutta.
 */

/**
 * Runge-Kutta integration step.
 */
template<typename t_simulation>
void algorithm::runge_kutta<t_simulation>
::step(t_simulation &sim,
       interaction::Forcefield<t_simulation> &ff,
       double const dt)
{
  // we do not have something as the old forces
  // so we do not exchange forces...

  int size = sim.system().pos().size();
  dx1.resize(size);
  dx2.resize(size);
  dx3.resize(size);
  dx4.resize(size);
  dv1.resize(size);
  dv2.resize(size);
  dv3.resize(size);
  dv4.resize(size);
  
  sim.system().exchange_vel();
  
  // first part
  dx1 = sim.system().old_vel();
  ff.calculate_interactions(sim);
  dv1 = sim.system().force() / sim.topology().mass();
  
  // second part
  sim.system().exchange_pos();
  dx2 = sim.system().old_vel() + dt/2.0 * dv1;
  sim.system().pos() = sim.system().old_pos() + dt/2.0 * dx1;
  ff.calculate_interactions(sim);
  dv2 = sim.system().force() / sim.topology().mass();
  
  // third part
  dx3 = sim.system().old_vel() + dt/2.0 * dv2;
  sim.system().pos() = sim.system().old_pos() + dt/2.0 * dx2;
  ff.calculate_interactions(sim);
  dv3 = sim.system().force() / sim.topology().mass();
  
  // fourth part
  dx4 = sim.system().old_vel() + dt * dv3;
  sim.system().pos() = sim.system().old_pos() + dt * dx3;
  ff.calculate_interactions(sim);
  dv4 = sim.system().force() / sim.topology().mass();
    
  // get the new positions / velocities by averageing (simpsons rule)
  sim.system().pos() = sim.system().old_pos() +
    dt/6.0 * (dx1 + 2*(dx2 + dx3) + dx4);
  
  sim.system().vel() = sim.system().old_vel() +
    dt/6.0 * (dv1 + 2*(dv2 + dv3) + dv4);
  
  // that's it
}
