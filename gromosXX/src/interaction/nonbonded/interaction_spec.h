/**
 * @file interaction_spec.h
 * specification of how to calculate the interactions
 */

#ifndef INCLUDED_INTERACTION_SPEC_H
#define INCLUDED_INTERACTION_SPEC_H

namespace interaction
{
  /**
   * @class Interaction_Spec
   * interaction specifications.
   */
  template<
    math::boundary_enum t_boundary = math::rectangular,
    bool t_perturbation = false,
    math::virial_enum t_virial = math::molecular_virial,
    bool t_atomic_cutoff = false,
    bool t_scaling = false
    >
  class Interaction_Spec
  {
  public:
    typedef Interaction_Spec<t_boundary,
			     t_perturbation,
			     t_virial, 
			     t_atomic_cutoff,
			     t_scaling>
    interaction_spec_type;

    static const math::boundary_enum boundary_type = t_boundary;
    static const bool do_exclusion = true;
    static const bool do_perturbation = t_perturbation;
    static const math::virial_enum do_virial = t_virial;
    static const bool do_atomic_cutoff = t_atomic_cutoff;
    static const bool do_scaling = t_scaling;
    static const bool do_bekker = false;
    
    // nonbonded interactions
    //==================================================
    typedef interaction::Exclusion_Filter<interaction_spec_type>
    exclusion_filter_type;
    
    typedef interaction::Range_Filter<interaction_spec_type>
    range_filter_type;
    
    typedef interaction::Nonbonded_Innerloop<interaction_spec_type>
    nonbonded_innerloop_type;
    
    typedef interaction::Standard_Pairlist_Algorithm<interaction_spec_type>
    pairlist_algorithm_type;
    
    typedef interaction::Nonbonded_Interaction<interaction_spec_type>
    nonbonded_interaction_type;

    // perturbed nonbonded interactions
    //==================================================

    typedef interaction::Perturbed_Nonbonded_Interaction<interaction_spec_type>
    perturbed_nonbonded_interaction_type;

    typedef interaction::Perturbed_Nonbonded_Innerloop<interaction_spec_type>
    perturbed_nonbonded_innerloop_type;

    typedef interaction::Perturbation_Filter<interaction_spec_type>
    perturbation_filter_type;
    
  };

 /**
   * @class Grid_Interaction_Spec
   * defines parameters for the interactions
   * using a grid based pairlist.
   */
  template<
    math::boundary_enum t_boundary = math::rectangular,
    bool t_perturbation = false,
    math::virial_enum t_virial = math::molecular_virial,
    bool t_atomic_cutoff = false,
    bool t_scaling = false>
  class Grid_Interaction_Spec
  {
  public:
    typedef Grid_Interaction_Spec<t_boundary,
				  t_perturbation,
				  t_virial,
				  t_atomic_cutoff,
				  t_scaling>
    interaction_spec_type;

    static const math::boundary_enum boundary_type = t_boundary;
    static const bool do_exclusion = true;
    static const bool do_perturbation = t_perturbation;
    static const math::virial_enum do_virial = t_virial;
    static const bool do_atomic_cutoff = t_atomic_cutoff;
    static const bool do_scaling = t_scaling;
    static const bool do_bekker = true;

    // nonbonded interactions
    //==================================================
    typedef interaction::Exclusion_Filter<interaction_spec_type>
    exclusion_filter_type;

    typedef interaction::Range_Filter<interaction_spec_type>
    range_filter_type;

    typedef interaction::Nonbonded_Innerloop<interaction_spec_type>
    nonbonded_innerloop_type;


    typedef interaction::Grid_Pairlist_Algorithm<interaction_spec_type>
    pairlist_algorithm_type;

    typedef interaction::Nonbonded_Interaction<interaction_spec_type>
    nonbonded_interaction_type;

    // perturbed nonbonded interactions
    //==================================================

    typedef interaction::Perturbed_Nonbonded_Interaction<interaction_spec_type>
    perturbed_nonbonded_interaction_type;

    typedef interaction::Perturbed_Nonbonded_Innerloop<interaction_spec_type>
    perturbed_nonbonded_innerloop_type;

    typedef interaction::Perturbation_Filter<interaction_spec_type>
    perturbation_filter_type;

  };
    
} // interaction

#endif
