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

    // and the interactions
    typedef interaction::Quartic_bond_interaction<simulation_type>
    quartic_bond_interaction_type;
    typedef interaction::harmonic_bond_interaction<simulation_type>
    harmonic_bond_interaction_type;
    typedef interaction::angle_interaction<simulation_type>
    angle_interaction_type;
    typedef interaction::Improper_dihedral_interaction<simulation_type>
    improper_interaction_type;
    typedef interaction::Dihedral_interaction<simulation_type>
    dihedral_interaction_type;
    
    // the nonbonded    
    //  - molecular virial
    //  - perturbation no
    //  - atomic cutoff no
    typedef interaction::Nonbonded_Interaction<
      simulation_type, 
      interaction::Nonbonded_Spec<simulation_type, 
				  interaction::molecular_virial,
				  false, false> >
    nonbonded_virial_interaction_type;

    //  - virial no
    //  - perturbation no
    //  - atomic cutoff no
    typedef interaction::Nonbonded_Interaction<
      simulation_type, 
      interaction::Nonbonded_Spec<simulation_type, 
				  interaction::no_virial,
				  false, false> >
    nonbonded_interaction_type;
    
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
    

    // and the interactions
    typedef interaction::Quartic_bond_interaction<simulation_type>
    quartic_bond_interaction_type;
    typedef interaction::harmonic_bond_interaction<simulation_type>
    harmonic_bond_interaction_type;
    typedef interaction::angle_interaction<simulation_type>
    angle_interaction_type;
    typedef interaction::Improper_dihedral_interaction<simulation_type>
    improper_interaction_type;
    typedef interaction::Dihedral_interaction<simulation_type>
    dihedral_interaction_type;
    
    // and the perturbed interactions
    typedef interaction::Perturbed_Quartic_Bond_Interaction<simulation_type>
    perturbed_quartic_bond_interaction_type;
    typedef interaction::Perturbed_Harmonic_Bond_Interaction<simulation_type>
    perturbed_harmonic_bond_interaction_type;
    typedef interaction::Perturbed_Angle_Interaction<simulation_type>
    perturbed_angle_interaction_type;
    typedef interaction::Perturbed_Improper_Dihedral_Interaction<simulation_type>
    perturbed_improper_interaction_type;
    typedef interaction::Perturbed_Dihedral_Interaction<simulation_type>
    perturbed_dihedral_interaction_type;
    
    // the nonbonded    
    //  - molecular virial
    //  - perturbation yes
    //  - atomic cutoff no
    typedef interaction::Perturbed_Nonbonded_Interaction<
      simulation_type, interaction::Nonbonded_Spec<simulation_type, interaction::molecular_virial,
						   true, false> >
    nonbonded_virial_interaction_type;

    typedef interaction::Perturbed_Nonbonded_Interaction<
      simulation_type, interaction::Nonbonded_Spec<simulation_type, interaction::no_virial,
						   true, false> >
    nonbonded_interaction_type;

  };

} // program


#endif
