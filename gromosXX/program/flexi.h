/**
 * @file flexi.h
 * specification of the typenames for various classes.
 */

#ifndef INCLUDED_FLEXI_H
#define INCLUDED_FLEXI_H

namespace program
{

  /**
   * @class FLEXI_spec
   * typedef's for the various classes
   * needed for an MD simulation.
   */
  class FLEXI_spec
  {
    // the standard types
  public:
    // the simulation
    typedef simulation::Simulation<simulation::Topology,
				   simulation::System<math::any> >
    simulation_type;
    
    typedef algorithm::Berendsen_Thermostat   temperature_type;

    typedef algorithm::Berendsen_Barostat     pressure_type;

    typedef algorithm::Flexible_Constraint<simulation_type> 
    distance_constraint_type;

    typedef algorithm::Leap_Frog<simulation_type>
    integration_type;    

  };

  /**
   * @class perturbed_FLEXI_spec
   * typedef's for the various classes
   * needed for an MD simulation with perturbation.
   */
  class perturbed_FLEXI_spec
  {
    // the standard types
  public:
    // the simulation
    typedef simulation::Simulation<simulation::Perturbation_Topology,
				   simulation::System<math::any> >
    simulation_type;
    
    typedef algorithm::Berendsen_Thermostat   temperature_type;

    typedef algorithm::Berendsen_Barostat     pressure_type;

    typedef algorithm::Perturbed_Flexible_Constraint<simulation_type> 
    distance_constraint_type;

    typedef algorithm::Leap_Frog<simulation_type>
    integration_type;
    
  };

} // program


#endif
