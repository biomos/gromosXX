/**
 * @file flexi.h
 * specification of the typenames for various classes.
 */

#ifndef INCLUDED_FASTI_H
#define INCLUDED_FASTI_H

namespace program
{

  /**
   * @class Fast_Nonbonded_Spec
   * defines parameters for the nonbonded interactions.
   */
  template<typename t_simulation, 
	   interaction::virial_enum t_virial = interaction::molecular_virial,
	   bool t_perturbation = false, bool t_atomic_cutoff = false>
  class Fast_Nonbonded_Spec
  {
  public:
    static const bool do_exclusion = true;
    static const bool do_perturbation = t_perturbation;
    static const interaction::virial_enum do_virial = t_virial;
    static const bool do_atomic_cutoff = t_atomic_cutoff;
    static const bool do_scaling = false;

    // the simulation
    typedef t_simulation simulation_type;
    // the spec
    typedef Fast_Nonbonded_Spec<simulation_type, t_virial,
				t_perturbation, t_atomic_cutoff> 
    nonbonded_spec_type;
    
    // exclusion filter
    typedef interaction::Exclusion_Filter<simulation_type, 
					  nonbonded_spec_type>
    exclusion_filter_type;

    // perturbation filter
    typedef interaction::Perturbation_Filter<simulation_type, 
					     nonbonded_spec_type>
    perturbation_filter_type;

    // range filter
    typedef interaction::Range_Filter<simulation_type,
				      nonbonded_spec_type>
    range_filter_type;
    
    // nonbonded innerloop (non perturbed interactions)
    typedef interaction::Nonbonded_Innerloop<simulation_type,
					     nonbonded_spec_type>
    nonbonded_innerloop_type;

    // nonbonded innnerloop (perturbed interactions)
    typedef interaction::Perturbed_Nonbonded_Innerloop<simulation_type,
						       nonbonded_spec_type>
    perturbed_nonbonded_innerloop_type;

    // pairlist algorithm (should know the correct
    // Nonbonded_Interaction or always assume
    // Perturbed_Nonbonded_Interaction??)
    typedef interaction::Grid_Pairlist_Algorithm<simulation_type,
						 nonbonded_spec_type>
    pairlist_algorithm_type;
    
  };


  /**
   * @class Fasti_spec
   * typedef's for the various classes
   * needed for an MD simulation.
   */
  class Fasti_spec
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
      Fast_Nonbonded_Spec<simulation_type, 
			  interaction::molecular_virial,
			  false, false> >
    nonbonded_virial_interaction_type;

    //  - virial no
    //  - perturbation no
    //  - atomic cutoff no
    typedef interaction::Nonbonded_Interaction<
      simulation_type, 
      Fast_Nonbonded_Spec<simulation_type, 
			  interaction::no_virial,
			  false, false> >
    nonbonded_interaction_type;
    
  };

  /**
   * @class perturbed_Fasti_spec
   * typedef's for the various classes
   * needed for an MD simulation with perturbation.
   */
  class perturbed_Fasti_spec
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
      simulation_type, Fast_Nonbonded_Spec<simulation_type, interaction::molecular_virial,
						   true, false> >
    nonbonded_virial_interaction_type;

    typedef interaction::Perturbed_Nonbonded_Interaction<
      simulation_type, Fast_Nonbonded_Spec<simulation_type, interaction::no_virial,
					   true, false> >
    nonbonded_interaction_type;

  };

} // program


#endif
