/**
 * @file nonbonded_spec.h
 * the nonbonded spec file.
 */

#ifndef INCLUDED_NONBONDED_SPEC_H
#define INCLUDED_NONBONDED_SPEC_H

namespace interaction
{
  /**
   * @class Nonbonded_Spec
   * defines parameters for the nonbonded interactions.
   */
  template<typename t_simulation, virial_enum t_virial = molecular_virial,
	   bool t_perturbation = false, bool t_atomic_cutoff = false>
  class Nonbonded_Spec
  {
  public:    
    static const bool do_exclusion = true;
    static const bool do_perturbation = t_perturbation;
    static const virial_enum do_virial = t_virial;
    static const bool do_atomic_cutoff = t_atomic_cutoff;
    static const bool do_scaling = false;

    // the simulation
    typedef t_simulation simulation_type;
    // the spec
    typedef Nonbonded_Spec<simulation_type, t_virial,
			   t_perturbation, t_atomic_cutoff> 
    nonbonded_spec_type;
    
    // exclusion filter
    typedef Exclusion_Filter<simulation_type, 
			     nonbonded_spec_type>
    exclusion_filter_type;

    // perturbation filter
    typedef Perturbation_Filter<simulation_type, 
				nonbonded_spec_type>
    perturbation_filter_type;

    // range filter
    typedef Range_Filter<simulation_type,
			 nonbonded_spec_type>
    range_filter_type;
    
    // nonbonded innerloop (non perturbed interactions)
    typedef Nonbonded_Innerloop<simulation_type,
				nonbonded_spec_type>
    nonbonded_innerloop_type;

    // nonbonded innnerloop (perturbed interactions)
    typedef Perturbed_Nonbonded_Innerloop<simulation_type,
					  nonbonded_spec_type>
    perturbed_nonbonded_innerloop_type;

    // pairlist algorithm (should know the correct
    // Nonbonded_Interaction or always assume
    // Perturbed_Nonbonded_Interaction??)
    typedef Standard_Pairlist_Algorithm<simulation_type,
					nonbonded_spec_type>
    pairlist_algorithm_type;
    
  };

  /**
   * @class Nonbonded_Scaled_Spec
   * defines parameters for the 
   * nonbonded interactions (with
   * scaling of interactions)
   */
  template<typename t_simulation, virial_enum t_virial = molecular_virial,
	   bool t_perturbation = false, bool t_atomic_cutoff = false>
  class Nonbonded_Scaled_Spec
  {
  public:    
    static const bool do_exclusion = true;
    static const bool do_perturbation = t_perturbation;
    static const virial_enum do_virial = t_virial;
    static const bool do_atomic_cutoff = t_atomic_cutoff;
    static const bool do_scaling = true;

    // the simulation
    typedef t_simulation simulation_type;
    // the spec
    typedef Nonbonded_Scaled_Spec<simulation_type, t_virial,
				  t_perturbation, t_atomic_cutoff> 
    nonbonded_spec_type;
    
    // exclusion filter
    typedef Exclusion_Filter<simulation_type, 
			     nonbonded_spec_type>
    exclusion_filter_type;

    // perturbation filter
    typedef Perturbation_Filter<simulation_type, 
				nonbonded_spec_type>
    perturbation_filter_type;

    // range filter
    typedef Range_Filter<simulation_type,
			 nonbonded_spec_type>
    range_filter_type;
    
    // nonbonded innerloop (non perturbed interactions)
    typedef Nonbonded_Innerloop<simulation_type,
				nonbonded_spec_type>
    nonbonded_innerloop_type;

    // nonbonded innnerloop (perturbed interactions)
    typedef Perturbed_Nonbonded_Innerloop<simulation_type,
					  nonbonded_spec_type>
    perturbed_nonbonded_innerloop_type;

    // pairlist algorithm (should know the correct
    // Nonbonded_Interaction or always assume
    // Perturbed_Nonbonded_Interaction??)
    typedef Standard_Pairlist_Algorithm<simulation_type,
					nonbonded_spec_type>
    pairlist_algorithm_type;
    
  };
  
} // interaction

#endif
