/**
 * @file innerloop_template.h
 * call a (member) function with the correct template parameters to
 * define an Nonbonded_Innerloop or Perturbed_Nonbonded_Innerloop class
 * with the correct boundary conditions, virial and interaction term
 * enumeration values.
 */

/**
 * call a function f with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function to construct a
 * Nonbonded_Innerloop : split the interaction term function
 */
#define SPLIT_INTERACTION_FUNC(f, bound, ...) \
  switch(sim.param().force.interaction_function){ \
    case simulation::lj_crf_func : \
      if (!sim.param().qmmm.dynamic_buffer_charges){ \
        f<Interaction_Spec<bound, simulation::lj_crf_func> >(__VA_ARGS__); \
      } else { \
        f<Interaction_Spec<bound, simulation::lj_crf_func, simulation::qm_buffer_charge> >(__VA_ARGS__); \
      } \
      break; \
    case simulation::lj_ls_func : \
      f<Interaction_Spec<bound, simulation::lj_ls_func> >(__VA_ARGS__); \
      break; \
    case simulation::lj_func : \
      f<Interaction_Spec<bound, simulation::lj_func> >(__VA_ARGS__); \
      break; \
    case simulation::cgrain_func : \
      f<Interaction_Spec<bound, simulation::cgrain_func> >(__VA_ARGS__); \
      break; \
    case simulation::cggromos_func : \
      f<Interaction_Spec<bound, simulation::cggromos_func> >(__VA_ARGS__); \
      break; \
    case simulation::pol_lj_crf_func : \
      if (sim.param().polarise.damp) { \
        switch(sim.param().polarise.efield_site) { \
          case simulation::ef_atom : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_lj_crf_func, pol_damping_on, simulation::ef_atom > >(__VA_ARGS__); \
            break; \
          case simulation::ef_cos : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_lj_crf_func, pol_damping_on, simulation::ef_cos > >(__VA_ARGS__); \
            break; \
          default: io::messages.add("Electric field calculation site not implemented.", "innerloop_template", io::message::error); \
        } \
      } else { \
        switch(sim.param().polarise.efield_site) { \
          case simulation::ef_atom : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_lj_crf_func, pol_damping_off, simulation::ef_atom > >(__VA_ARGS__); \
            break; \
          case simulation::ef_cos : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_lj_crf_func, pol_damping_off, simulation::ef_cos > >(__VA_ARGS__); \
            break; \
          default: io::messages.add("Electric field calculation site not implemented.", "innerloop_template", io::message::error); \
        } \
      } \
      break; \
    case simulation::pol_off_lj_crf_func : \
      if (sim.param().polarise.damp) { \
        switch(sim.param().polarise.efield_site) { \
          case simulation::ef_atom : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_off_lj_crf_func, pol_damping_on, simulation::ef_atom > >(__VA_ARGS__); \
            break; \
          case simulation::ef_cos : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_off_lj_crf_func, pol_damping_on, simulation::ef_cos > >(__VA_ARGS__); \
            break; \
          default: io::messages.add("Electric field calculation site not implemented.", "innerloop_template", io::message::error); \
        } \
      } else { \
        switch(sim.param().polarise.efield_site) { \
          case simulation::ef_atom : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_off_lj_crf_func, pol_damping_off, simulation::ef_atom > >(__VA_ARGS__); \
            break; \
          case simulation::ef_cos : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_off_lj_crf_func, pol_damping_off, simulation::ef_cos > >(__VA_ARGS__); \
            break; \
          default: io::messages.add("Electric field calculation site not implemented.", "innerloop_template", io::message::error); \
        } \
      } \
      break; \
    default: \
      io::messages.add("wrong interaction function", "innerloop_template", io::message::error); \
      \
  } \


/**
 * call a function with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function for a
 * Nonbonded_Innerloop : split the boundary
 *
 */
#define SPLIT_INNERLOOP(f, ...) \
  switch(conf.boundary_type){ \
    case math::vacuum : \
      SPLIT_INTERACTION_FUNC(f, math::vacuum, __VA_ARGS__); \
      break; \
    case math::rectangular : \
      SPLIT_INTERACTION_FUNC(f, math::rectangular, __VA_ARGS__); \
      break; \
    case math::truncoct : \
    case math::triclinic : \
      SPLIT_INTERACTION_FUNC(f, math::triclinic, __VA_ARGS__); \
      break; \
    default: \
      io::messages.add("wrong boundary type", "template_split", io::message::error); \
  } \

/**
 * call a function with specified interaction function using the correct values for
 * boundary, virial and interaction term function for a
 * Nonbonded_Innerloop : split the boundary
 *
 */
#define SPLIT_MY_INNERLOOP(f, interaction_function, ...) \
  switch(conf.boundary_type){ \
    case math::vacuum : \
      f<Interaction_Spec<math::vacuum, interaction_function> >(__VA_ARGS__); \
      break; \
    case math::rectangular : \
      f<Interaction_Spec<math::rectangular, interaction_function> >(__VA_ARGS__); \
      break; \
    case math::truncoct : \
    case math::triclinic : \
      f<Interaction_Spec<math::triclinic, interaction_function> >(__VA_ARGS__); \
      break; \
    default: \
      io::messages.add("wrong boundary type", "template_split", io::message::error); \
  } \
  
////////////////////////////////////////////////////////////////////////////////

/**
 * call a function f with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function to construct a
 * Nonbonded_Innerloop : split the interaction term function
 */
#define PERT_SPLIT_INTERACTION_FUNC(f, pertspec, bound, ...) \
  switch(sim.param().force.interaction_function){ \
    case simulation::lj_crf_func : \
      f< Interaction_Spec<bound, simulation::lj_crf_func>, \
         pertspec > (__VA_ARGS__); \
      break; \
    case simulation::cgrain_func : \
      f< Interaction_Spec<bound, simulation::cgrain_func>, \
         pertspec > (__VA_ARGS__); \
      break; \
    case simulation::cggromos_func : \
      f< Interaction_Spec<bound, simulation::cggromos_func>, \
         pertspec > (__VA_ARGS__); \
      break; \
    case simulation::pol_lj_crf_func : \
      if (sim.param().polarise.damp) { \
        switch(sim.param().polarise.efield_site) { \
          case simulation::ef_atom : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_lj_crf_func, pol_damping_on, simulation::ef_atom >, pertspec>(__VA_ARGS__); \
            break; \
          case simulation::ef_cos : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_lj_crf_func, pol_damping_on, simulation::ef_cos >, pertspec>(__VA_ARGS__); \
            break; \
          default: io::messages.add("Electric field calculation site not implemented.", "innerloop_template", io::message::error); \
        } \
      } else { \
        switch(sim.param().polarise.efield_site) { \
          case simulation::ef_atom : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_lj_crf_func, pol_damping_off, simulation::ef_atom >, pertspec>(__VA_ARGS__); \
            break; \
          case simulation::ef_cos : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_lj_crf_func, pol_damping_off, simulation::ef_cos >, pertspec>(__VA_ARGS__); \
            break; \
          default: io::messages.add("Electric field calculation site not implemented.", "innerloop_template", io::message::error); \
        } \
      } \
      break; \
    case simulation::pol_off_lj_crf_func : \
      if (sim.param().polarise.damp) { \
        switch(sim.param().polarise.efield_site) { \
          case simulation::ef_atom : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_off_lj_crf_func, pol_damping_on, simulation::ef_atom >, pertspec>(__VA_ARGS__); \
            break; \
          case simulation::ef_cos : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_off_lj_crf_func, pol_damping_on, simulation::ef_cos >, pertspec>(__VA_ARGS__); \
            break; \
          default: io::messages.add("Electric field calculation site not implemented.", "innerloop_template", io::message::error); \
        } \
      } else { \
        switch(sim.param().polarise.efield_site) { \
          case simulation::ef_atom : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_off_lj_crf_func, pol_damping_off, simulation::ef_atom >, pertspec>(__VA_ARGS__); \
            break; \
          case simulation::ef_cos : \
            f<Polarisation_Interaction_Spec<bound, simulation::pol_off_lj_crf_func, pol_damping_off, simulation::ef_cos >, pertspec>(__VA_ARGS__); \
            break; \
          default: io::messages.add("Electric field calculation site not implemented.", "innerloop_template", io::message::error); \
        } \
      } \
      break; \
    default: \
      io::messages.add("wrong interaction function", "innerloop_template", io::message::error); \
      \
  } \


/**
 * call a function with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function for a
 * Nonbonded_Innerloop : split the boundary
 */
#define PERT_SPLIT_BOUNDARY(f, pertspec, ...) \
  switch(conf.boundary_type){ \
    case math::vacuum : \
      PERT_SPLIT_INTERACTION_FUNC(f, pertspec, math::vacuum, __VA_ARGS__); \
      break; \
    case math::rectangular : \
      PERT_SPLIT_INTERACTION_FUNC(f, pertspec, math::rectangular, __VA_ARGS__); \
      break; \
    case math::truncoct : \
    case math::triclinic : \
      PERT_SPLIT_INTERACTION_FUNC(f, pertspec, math::triclinic, __VA_ARGS__); \
      break; \
    default: \
      io::messages.add("wrong boundary type", "template_split", io::message::error); \
  } \

/**
 * call a function with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function for a
 * Nonbonded_Innerloop : split perturbation
 */
#define SPLIT_PERT_INNERLOOP(f, ...) \
  assert(sim.param().perturbation.perturbation || sim.param().eds.eds); \
  if (sim.param().perturbation.scaling){ \
    PERT_SPLIT_BOUNDARY(f, Perturbation_Spec<scaling_on>, __VA_ARGS__); \
    } \
  else { \
    PERT_SPLIT_BOUNDARY(f, Perturbation_Spec<scaling_off>, __VA_ARGS__); \
  } \


////////////////////////////////////////////////////////////////////////////////

/**
 * call a function with the appropriate values for scaling
 * (assuming that perturbation is enabled).
 */
#define SPLIT_PERTURBATION(f, ...) \
assert(sim.param().perturbation.perturbation); \
if (sim.param().perturbation.scaling){ \
  f<Perturbation_Spec<scaling_on> >(__VA_ARGS__); \
  } else { \
  f<Perturbation_Spec<scaling_off> >(__VA_ARGS__); \
} \

////////////////////////////////////////////////////////////////////////////////

