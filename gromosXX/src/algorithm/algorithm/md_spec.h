/**
 * @file md_spec.h
 * specification of the typenames for various classes.
 */

#ifndef INCLUDED_MD_SPEC_H
#define INCLUDED_MD_SPEC_H

namespace algorithm
{

  /**
   * @class MD_spec
   * typedef's for the various classes
   * needed for an MD simulation.
   */
  class MD_spec
  {
    // the standard types
  public:
    // the simulation
    typedef simulation::Simulation<simulation::Topology,
				   simulation::System<math::any> >
    simulation_type;
    
    typedef algorithm::Berendsen_Thermostat   temperature_type;
    typedef algorithm::Berendsen_Barostat     pressure_type;
    typedef algorithm::Shake<simulation_type> 
    distance_constraint_type;
    typedef algorithm::Leap_Frog<simulation_type>
    integration_type;
  };

  /**
   * @class perturbed_MD_spec
   * typedef's for the various classes
   * needed for an MD simulation with perturbation.
   */
  class perturbed_MD_spec
  {
    // the standard types
  public:
    // the simulation
    typedef simulation::Simulation<simulation::Perturbation_Topology,
				   simulation::System<math::any> >
    simulation_type;
    
    typedef algorithm::Berendsen_Thermostat   temperature_type;
    typedef algorithm::Berendsen_Barostat     pressure_type;
    typedef algorithm::Perturbed_Shake<simulation_type> 
    distance_constraint_type;
    typedef algorithm::Leap_Frog<simulation_type>
    integration_type;
    
  };
  
} // algorithm


#endif
