/**
 * @file innerloop_template.h
 * call a (member) function with the correct template parameters to
 * define an Nonbonded_Innerloop or Perturbed_Nonbonded_Innerloop class
 * with the correct boundary conditions, virial and interaction term
 * enumeration values.
 */

/**
 * @define SPLIT_INTERACTION_FUNC
 * call a function f with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function to construct a
 * Nonbonded_Innerloop : split the interaction term function
 */
#define SPLIT_INTERACTION_FUNC(f, bound, vir, ...) \
  switch(sim.param().force.interaction_function){ \
    case simulation::lj_crf_func : \
      f<Interaction_Spec<bound, vir, simulation::lj_crf_func> >(__VA_ARGS__); \
      break; \
    default: \
      io::messages.add("wrong interaction function", "innerloop_template", io::message::error); \
      \
  } \

/**
 * @define SPLIT_VIRIAL
 * call a function f with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function to construct a
 * Nonbonded_Innerloop : split the virial
 */
#define SPLIT_VIRIAL(f, bound, ...) \
  switch(sim.param().pcouple.virial){ \
    case math::no_virial : \
      SPLIT_INTERACTION_FUNC(f, bound, math::no_virial, __VA_ARGS__); \
      break; \
    case math::molecular_virial : \
      SPLIT_INTERACTION_FUNC(f, bound, math::molecular_virial, __VA_ARGS__); \
      break; \
    case math::atomic_virial : \
      SPLIT_INTERACTION_FUNC(f, bound, math::atomic_virial, __VA_ARGS__); \
      break; \
    default: \
      io::messages.add("wrong virial type", "template_split", io::message::error); \
      \
  } \


/**
 * @define SPLIT_INNERLOOP
 * call a function with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function for a
 * Nonbonded_Innerloop : split the boundary
 *
 */
#define SPLIT_INNERLOOP(f, ...) \
  switch(conf.boundary_type){ \
    case math::vacuum : \
      SPLIT_INTERACTION_FUNC(f, math::vacuum, math::no_virial, __VA_ARGS__); \
      break; \
    case math::rectangular : \
      SPLIT_VIRIAL(f, math::rectangular, __VA_ARGS__); \
      break; \
    case math::truncoct : \
      SPLIT_VIRIAL(f, math::truncoct, __VA_ARGS__); \
      break; \
    default: \
      io::messages.add("wrong boundary type", "template_split", io::message::error); \
  } \

////////////////////////////////////////////////////////////////////////////////

/**
 * @define PERT_SPLIT_INTERACTION_FUNC
 * call a function f with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function to construct a
 * Nonbonded_Innerloop : split the interaction term function
 */
#define PERT_SPLIT_INTERACTION_FUNC(f, bound, vir, ...) \
  switch(sim.param().force.interaction_function){ \
    case simulation::lj_crf_func : \
      f<Interaction_Spec<bound, vir, simulation::lj_crf_func>, \
        t_perturbation_details> (__VA_ARGS__); break; \
      break; \
    default: \
      io::messages.add("wrong interaction function", "innerloop_template", io::message::error); \
      \
  } \


/**
 * @define PERT_SPLIT_VIRIAL
 * call a function f with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function to construct a
 * Nonbonded_Innerloop : split the virial
 */
#define PERT_SPLIT_VIRIAL(f, bound, ...) \
  switch(sim.param().pcouple.virial){ \
    case math::no_virial : \
      PERT_SPLIT_INTERACTION_FUNC(f, bound, math::no_virial, __VA_ARGS__); \
      break; \
    case math::molecular_virial : \
      PERT_SPLIT_INTERACTION_FUNC(f, bound, math::molecular_virial, __VA_ARGS__); \
      break; \
    case math::atomic_virial : \
      PERT_SPLIT_INTERACTION_FUNC(f, bound, math::atomic_virial, __VA_ARGS__); \
      break; \
    default: \
      io::messages.add("wrong virial type", "template_split", io::message::error); \
      \
  } \


/**
 * @define PERT_SPLIT_PERT_INNERLOOP
 * call a function with a Interaction_Spec using the correct values for
 * boundary, virial and interaction term function for a
 * nonbonded innerloop.
 * the function also gets a second template typename, t_perturbation_details.
 * that has to be defined already (you should use SPLIT_PERTURBATION before).
 */
#define SPLIT_PERT_INNERLOOP(f, ...) \
  switch(conf.boundary_type){ \
    case math::vacuum : \
      PERT_SPLIT_INTERACTION_FUNC(f, math::vacuum, math::no_virial, __VA_ARGS__); \
      break; \
    case math::rectangular : \
      PERT_SPLIT_VIRIAL(f, math::rectangular, __VA_ARGS__); \
      break; \
    case math::truncoct : \
      PERT_SPLIT_VIRIAL(f, math::truncoct, __VA_ARGS__); \
      break; \
    default: \
      io::messages.add("wrong boundary type", "template_split", io::message::error); \
  } \


////////////////////////////////////////////////////////////////////////////////

/**
 * @define SPLIT_PERTURBATION
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

