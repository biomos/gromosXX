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
    typedef interaction::Basic_Pairlist<
      simulation_type, interaction::Chargegroup_Range_Pairlist_Algorithm<
      simulation_type, interaction::Twinrange_Chargegroup_Filter<
      simulation_type, interaction::Nonbonded_Base, 
      interaction::Nonbonded_Inner_Loop<
      simulation_type, interaction::Storage> 
    > > > 
    pairlist_type;
      
    typedef interaction::Nonbonded_Inner_Loop<
      simulation_type, simulation_type::system_type>
    innerloop_type;

    typedef interaction::Basic_Pairlist<
      simulation_type, interaction::Chargegroup_Range_Pairlist_Algorithm<
      simulation_type, interaction::Twinrange_Chargegroup_Filter<
      simulation_type, interaction::Nonbonded_Base, 
      interaction::Nonbonded_Inner_Loop_Virial<
      simulation_type, interaction::Storage
      > > > > 
    pairlist_virial_type;
      
    typedef interaction::Nonbonded_Inner_Loop_Virial<
      simulation_type, simulation_type::system_type> innerloop_virial_type;
    
    typedef interaction::Nonbonded_Interaction<
      simulation_type, pairlist_virial_type, innerloop_virial_type>
    nonbonded_virial_interaction_type;

    typedef interaction::Nonbonded_Interaction<
      simulation_type, pairlist_type, innerloop_type>
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
    
    // the nonbonded    
    // both nonbonded use the same (perturbed) pairlist !!!
    // and inner loop!!!
    typedef interaction::Perturbed_Nonbonded_Inner_Loop<
      simulation_type, simulation_type::system_type>
    innerloop_type;

    typedef interaction::Perturbed_Nonbonded_Inner_Loop_Virial<
      simulation_type, simulation_type::system_type>
    innerloop_virial_type;     

    typedef interaction::Basic_Pairlist<
      simulation_type, interaction::Chargegroup_Range_Pairlist_Algorithm<
      simulation_type, interaction::Twinrange_Chargegroup_Filter<
      simulation_type, interaction::Nonbonded_Base, 
      interaction::Perturbed_Nonbonded_Inner_Loop<
      simulation_type, interaction::Storage>, 
      interaction::Perturbation_Filter<
      simulation_type,interaction::Nonbonded_Base, true
      > > > >
    pairlist_type;

    typedef interaction::Basic_Pairlist<
      simulation_type, interaction::Chargegroup_Range_Pairlist_Algorithm<
      simulation_type, interaction::Twinrange_Chargegroup_Filter<
      simulation_type, interaction::Nonbonded_Base,
      interaction::Perturbed_Nonbonded_Inner_Loop_Virial<
      simulation_type, interaction::Storage>,
      interaction::Perturbation_Filter<
      simulation_type, interaction::Nonbonded_Base, true
      > > > > 
    pairlist_virial_type;

    typedef interaction::Nonbonded_Interaction<simulation_type,
					       pairlist_virial_type,
					       innerloop_virial_type>
    nonbonded_virial_interaction_type;

    typedef interaction::Nonbonded_Interaction<simulation_type,
					       pairlist_type,
					       innerloop_type>
    nonbonded_interaction_type;
       

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
    
    // the perturbed nonbonded
    typedef interaction::Perturbed_Nonbonded_Interaction<
      simulation_type, pairlist_virial_type, innerloop_virial_type,
      nonbonded_virial_interaction_type
      >
    perturbed_nonbonded_virial_interaction_type;

    typedef interaction::Perturbed_Nonbonded_Interaction<
      simulation_type, pairlist_type, innerloop_type,
      nonbonded_interaction_type
      >
    perturbed_nonbonded_interaction_type; 
  };

} // program


#endif
