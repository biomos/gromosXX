/**
 * @file template_split.h
 * split function call based on template paramters
 */

#define SPLIT_BOUNDARY(f, ...) \
switch(conf.boundary_type){ \
  case math::vacuum : f<math::vacuum>(__VA_ARGS__); break; \
  case math::rectangular : f<math::rectangular>(__VA_ARGS__); break; \
  case math::triclinic : f<math::triclinic>(__VA_ARGS__); break; \
  default: throw std::string("wrong boundary type"); \
} \


#define SPLIT_VIRIAL_BOUNDARY(f, ...) \
switch(conf.boundary_type){ \
  case math::vacuum : f<math::vacuum, math::no_virial>(__VA_ARGS__); break; \
  case math::rectangular : \
    switch(sim.param().pcouple.virial){ \
      case math::no_virial : f<math::rectangular, math::no_virial>(__VA_ARGS__); break; \
      case math::molecular_virial : f<math::rectangular, math::molecular_virial>(__VA_ARGS__); break; \
      case math::atomic_virial : f<math::rectangular, math::atomic_virial>(__VA_ARGS__); break; \
      default: throw std::string("wrong virial type"); \
    } \
    break; \
  default: throw std::string("wrong boundary type"); \
} \


#define SPLIT_PERTURBATION(f, ...) \
assert(sim.param().perturbation.perturbation); \
if (sim.param().perturbation.scaling){ \
  f<Perturbation_Spec<scaling_on> >(__VA_ARGS__); \
} else { \
  f<Perturbation_Spec<scaling_off> >(__VA_ARGS__); \
} \


#define SPLIT_INNERLOOP(f, ...) \
if (sim.param().pairlist.grid) { \
  switch(conf.boundary_type){ \
    case math::vacuum : f<Interaction_Spec<math::vacuum, math::no_virial, true> >(__VA_ARGS__); break; \
    case math::rectangular : \
      switch(sim.param().pcouple.virial){ \
        case math::no_virial : f<Interaction_Spec<math::rectangular, math::no_virial, true> >(__VA_ARGS__); break; \
        case math::molecular_virial : f<Interaction_Spec<math::rectangular, math::molecular_virial, true> >(__VA_ARGS__); break; \
        case math::atomic_virial : f<Interaction_Spec<math::rectangular, math::atomic_virial, true> >(__VA_ARGS__); break; \
        default: throw std::string("wrong virial type"); \
      } \
      break; \
    default: throw std::string("wrong boundary type"); \
  } \
} \
else { \
  switch(conf.boundary_type){ \
    case math::vacuum : f<Interaction_Spec<math::vacuum, math::no_virial, false> >(__VA_ARGS__); break; \
    case math::rectangular : \
      switch(sim.param().pcouple.virial){ \
        case math::no_virial : f<Interaction_Spec<math::rectangular, math::no_virial, false> >(__VA_ARGS__); break; \
        case math::molecular_virial : f<Interaction_Spec<math::rectangular, math::molecular_virial, false> >(__VA_ARGS__); break; \
        case math::atomic_virial : f<Interaction_Spec<math::rectangular, math::atomic_virial, false> >(__VA_ARGS__); break; \
        default: throw std::string("wrong virial type"); \
      } \
      break; \
    default: throw std::string("wrong boundary type"); \
  } \
} \

#define SPLIT_INNERLOOP_NO_GRID(f, ...) \
switch(conf.boundary_type){ \
  case math::vacuum : f<Interaction_Spec<math::vacuum, math::no_virial, false> >(__VA_ARGS__); break; \
  case math::rectangular : \
    switch(sim.param().pcouple.virial){ \
      case math::no_virial : f<Interaction_Spec<math::rectangular, math::no_virial, false> >(__VA_ARGS__); break; \
      case math::molecular_virial : f<Interaction_Spec<math::rectangular, math::molecular_virial, false> >(__VA_ARGS__); break; \
      case math::atomic_virial : f<Interaction_Spec<math::rectangular, math::atomic_virial, false> >(__VA_ARGS__); break; \
      default: throw std::string("wrong virial type"); \
    } \
    break; \
  default: throw std::string("wrong boundary type"); \
} \


#define SPLIT_PERT_INNERLOOP(f, ...) \
if (sim.param().pairlist.grid) { \
  switch(conf.boundary_type){ \
    case math::vacuum : f<Interaction_Spec<math::vacuum, math::no_virial, true>, t_perturbation_details>(__VA_ARGS__); break; \
    case math::rectangular : \
      switch(sim.param().pcouple.virial){ \
        case math::no_virial : f<Interaction_Spec<math::rectangular, math::no_virial, true>, t_perturbation_details>(__VA_ARGS__); break; \
        case math::molecular_virial : f<Interaction_Spec<math::rectangular, math::molecular_virial, true>, t_perturbation_details>(__VA_ARGS__); break; \
        case math::atomic_virial : f<Interaction_Spec<math::rectangular, math::atomic_virial, true>, t_perturbation_details>(__VA_ARGS__); break; \
        default: throw std::string("wrong virial type"); \
      } \
      break; \
    default: throw std::string("wrong boundary type"); \
  } \
} \
else { \
  switch(conf.boundary_type){ \
    case math::vacuum : f<Interaction_Spec<math::vacuum, math::no_virial, false>, t_perturbation_details>(__VA_ARGS__); break; \
    case math::rectangular : \
      switch(sim.param().pcouple.virial){ \
        case math::no_virial : f<Interaction_Spec<math::rectangular, math::no_virial, false>, t_perturbation_details>(__VA_ARGS__); break; \
        case math::molecular_virial : f<Interaction_Spec<math::rectangular, math::molecular_virial, false>, t_perturbation_details>(__VA_ARGS__); break; \
        case math::atomic_virial : f<Interaction_Spec<math::rectangular, math::atomic_virial, false>, t_perturbation_details>(__VA_ARGS__); break; \
        default: throw std::string("wrong virial type"); \
      } \
      break; \
    default: throw std::string("wrong boundary type"); \
  } \
} \

