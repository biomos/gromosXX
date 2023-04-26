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
  const bool scaling_on = true;
  const bool scaling_off = false;
  const bool pol_damping_on = true;
  const bool pol_damping_off = false;
  
  /**
   * @class Interaction_Spec
   * interaction specifications.
   */
  template<
    math::boundary_enum t_boundary = math::rectangular,
    // math::virial_enum t_virial = math::molecular_virial,
    simulation::interaction_func_enum t_interaction_func = simulation::lj_crf_func,
    simulation::charge_type_enum t_charge_type = simulation::mm_charge,
    bool t_pol_damping = false,
    simulation::efield_site_enum t_efield_site = simulation::ef_atom
  >
  class Interaction_Spec
  {
  public:
    typedef Interaction_Spec<t_boundary,
			     // t_virial,
			     t_interaction_func
			     >
    interaction_spec_type;

    static const math::boundary_enum boundary_type = t_boundary;
    // static const math::virial_enum do_virial = t_virial;
    static const simulation::interaction_func_enum interaction_func = t_interaction_func;
    static const bool pol_damping = t_pol_damping;
    static const simulation::efield_site_enum efield_site = t_efield_site;
    static const simulation::charge_type_enum charge_type = t_charge_type;
  };

  /**
   * @class Polarisation_Interaction_Spec
   * adapter class for polarisation interaction specifications.
   */
  template<
    math::boundary_enum t_boundary = math::rectangular,
    simulation::interaction_func_enum t_interaction_func = simulation::lj_crf_func,
    bool t_pol_damping = false,
    simulation::efield_site_enum t_efield_site = simulation::ef_atom,
    simulation::charge_type_enum t_charge_type = simulation::mm_charge
  >
  class Polarisation_Interaction_Spec :
    public Interaction_Spec<t_boundary, t_interaction_func, t_charge_type, t_pol_damping, t_efield_site>
  {};

  /**
   * @class Perturbation_Spec
   * specifies if and what kind of perturbation to do (or rather apply)
   */
  template<
    bool t_scaling = scaling_off
  >
  class Perturbation_Spec
  {
  public:
    typedef Perturbation_Spec<t_scaling>
    perturbation_spec_type;
    
    static const bool do_scaling = t_scaling;
  };

} // interaction

#endif
