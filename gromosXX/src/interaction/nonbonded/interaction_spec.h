/**
 * @file interaction_spec.h
 * specification of how to calculate the interactions
 */

#ifndef INCLUDED_INTERACTION_SPEC_H
#define INCLUDED_INTERACTION_SPEC_H

namespace interaction
{
  const bool atomic_cutoff_on = true;
  const bool atomic_cutoff_off = false;
  const bool perturbation_on = true;
  const bool perturbation_off = false;
  const bool bekker_on = true;
  const bool bekker_off = false;
  const bool scaling_on = true;
  const bool scaling_off = false;
  
  /**
   * @class Interaction_Spec
   * interaction specifications.
   */
  template<
    math::boundary_enum t_boundary = math::rectangular,
    math::virial_enum t_virial = math::molecular_virial,
    bool t_atomic_cutoff = atomic_cutoff_off,
    bool t_bekker = bekker_off,
    bool t_scaling = scaling_off    
    >
  class Interaction_Spec
  {
  public:
    typedef Interaction_Spec<t_boundary,
			     t_virial, 
			     t_atomic_cutoff,
			     t_bekker,
			     t_scaling>
    interaction_spec_type;

    static const math::boundary_enum boundary_type = t_boundary;
    static const bool do_exclusion = true;
    static const math::virial_enum do_virial = t_virial;
    static const bool do_atomic_cutoff = t_atomic_cutoff;
    static const bool do_scaling = t_scaling;
    static const bool do_bekker = t_bekker;
    
  };

} // interaction

#endif
