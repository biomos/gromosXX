/**
 * @file interaction_spec.h
 * specification of how to calculate the interactions
 */

#ifndef INCLUDED_INTERACTION_SPEC_H
#define INCLUDED_INTERACTION_SPEC_H

namespace algorithm
{
  /**
   * @class Interaction_spec
   * interaction specifications.
   */
  template<typename t_simulation,
	   bool t_perturbation = false,
	   interaction::virial_enum t_virial 
	   = interaction::molecular_virial,
	   bool t_atomic_cutoff = false,
	   bool t_scaling = false>
  class Interaction_spec
  {
  public:
    typedef t_simulation simulation_type;
    typedef Interaction_spec<t_simulation, t_perturbation,
			     t_virial, t_atomic_cutoff,
			     t_scaling>
    interaction_spec_type;

    static const bool do_exclusion = true;
    static const bool do_perturbation = t_perturbation;
    static const interaction::virial_enum do_virial = t_virial;
    static const bool do_atomic_cutoff = t_atomic_cutoff;
    static const bool do_scaling = t_scaling;
    
    // the interaction types
    // bonded interactions
    //==============================
    typedef interaction::Quartic_bond_interaction<simulation_type, 
						  interaction_spec_type>
    quartic_bond_interaction_type;

    typedef interaction::harmonic_bond_interaction<simulation_type,
						   interaction_spec_type>
    harmonic_bond_interaction_type;

    typedef interaction::angle_interaction<simulation_type,
					   interaction_spec_type>
    angle_interaction_type;

    typedef interaction::Improper_dihedral_interaction<simulation_type,
						       interaction_spec_type>
    improper_interaction_type;

    typedef interaction::Dihedral_interaction<simulation_type,
					      interaction_spec_type>
    dihedral_interaction_type;
    
    // nonbonded interactions
    //==============================
    typedef interaction::Exclusion_Filter<simulation_type,
					  interaction_spec_type>
    exclusion_filter_type;
    
    typedef interaction::Perturbation_Filter<simulation_type,
					     interaction_spec_type>
    perturbation_filter_type;
    
    typedef interaction::Range_Filter<simulation_type,
				      interaction_spec_type>
    range_filter_type;
    
    typedef interaction::Nonbonded_Innerloop<simulation_type,
					     interaction_spec_type>
    nonbonded_innerloop_type;
    
    typedef interaction::Perturbed_Nonbonded_Innerloop<simulation_type,
						       interaction_spec_type>
    perturbed_nonbonded_innerloop_type;
    
    typedef interaction::Standard_Pairlist_Algorithm<simulation_type,
						     interaction_spec_type>
    pairlist_algorithm_type;
    
    typedef interaction::Nonbonded_Interaction<
      simulation_type, 
      interaction_spec_type
      >
    nonbonded_interaction_type;

    // perturbed interactions
    typedef interaction::Perturbed_Quartic_Bond_Interaction<simulation_type,
							    interaction_spec_type>
    perturbed_quartic_bond_interaction_type;

    typedef interaction::Perturbed_Harmonic_Bond_Interaction<simulation_type,
							     interaction_spec_type>
    perturbed_harmonic_bond_interaction_type;

    typedef interaction::Perturbed_Angle_Interaction<simulation_type,
						     interaction_spec_type>
    perturbed_angle_interaction_type;

    typedef interaction::
    Perturbed_Improper_Dihedral_Interaction<simulation_type,
					    interaction_spec_type>
    perturbed_improper_interaction_type;

    typedef interaction::Perturbed_Dihedral_Interaction<simulation_type,
							interaction_spec_type>
    perturbed_dihedral_interaction_type;

    // perturbed nonbonded interactions
    typedef interaction::Perturbed_Nonbonded_Interaction<
      simulation_type, 
      interaction_spec_type>
    perturbed_nonbonded_interaction_type;

  };
    
} // algorithm

#endif
