/**
 * @file template_split.h
 * split function call based on template paramters
 */

////////////////////////////////////////////////////////////////////////////////
/**
 * @define SPLIT_BOUNDARY
 * call a function with the appropriate boundary as template parameter
 */
#define SPLIT_BOUNDARY(f, ...) \
switch(conf.boundary_type){ \
  case math::vacuum : f<math::vacuum>(__VA_ARGS__); break; \
  case math::rectangular : f<math::rectangular>(__VA_ARGS__); break; \
  case math::truncoct : \
  case math::triclinic : f<math::triclinic>(__VA_ARGS__); break; \
  default: io::messages.add("wrong boundary type", "template_split", io::message::error); \
} \

////////////////////////////////////////////////////////////////////////////////
/**
 * @define SPLIT_MY_BOUNDARY
 * call a function with the appropriate boundary as template parameter
 */
#define SPLIT_MY_BOUNDARY(bound, f, ...) \
switch(bound){ \
  case math::vacuum : f<math::vacuum>(__VA_ARGS__); break; \
  case math::rectangular : f<math::rectangular>(__VA_ARGS__); break; \
  case math::truncoct : \
  case math::triclinic : f<math::triclinic>(__VA_ARGS__); break; \
  default: io::messages.add("wrong boundary type", "template_split", io::message::error); \
} \

////////////////////////////////////////////////////////////////////////////////

/**
 * @define SPLIT_VIRIAL_BOUNDARY
 * call a function with the appropriate values for virial and boundary
 * as template parameters.
 */
#define SPLIT_VIRIAL_BOUNDARY(f, ...) \
switch(conf.boundary_type){ \
  case math::vacuum : f<math::vacuum, math::no_virial>(__VA_ARGS__); break; \
  case math::rectangular : \
    switch(sim.param().pcouple.virial){ \
      case math::no_virial : f<math::rectangular, math::no_virial>(__VA_ARGS__); break; \
      case math::molecular_virial : f<math::rectangular, math::molecular_virial>(__VA_ARGS__); break; \
      case math::atomic_virial : f<math::rectangular, math::atomic_virial>(__VA_ARGS__); break; \
      default: io::messages.add("wrong virial type", "template_split", io::message::error); \
    } \
    break; \
  case math::truncoct : \
  case math::triclinic : \
    switch(sim.param().pcouple.virial){ \
      case math::no_virial : f<math::triclinic, math::no_virial>(__VA_ARGS__); break; \
      case math::molecular_virial : f<math::triclinic, math::molecular_virial>(__VA_ARGS__); break; \
      case math::atomic_virial : f<math::triclinic, math::atomic_virial>(__VA_ARGS__); break; \
      default: io::messages.add("wrong virial type", "template_split", io::message::error); \
    } \
    break; \
  default: io::messages.add("wrong boundary type", "template_split", io::message::error); \
} \

////////////////////////////////////////////////////////////////////////////////


