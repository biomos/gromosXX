/**
 * @file debug.h
 * define debug macros.
 */
 /**
 * @page debug Debug Levels
 *  -  0: no output
 *  -  1: important program decisions
 *  -  2: important program path
 *  -  3: 
 *  -  4: program path (not unimportant functions)
 *  -  5: function outline
 *  -  6:
 *  -  7: function detail
 *  -  8:
 *  -  9:
 *  - 10 : everything
 *  - 15 : even the ridiculous
 *
 * every file that uses it has to first
 * #define MODULE to the namespace name
 * #define SUBMODULE to the subdirectory in the namespace
 * #include debug.h (../../debug.h usually)
 * add the appropriate variables to namespace.h, namespace.cc
 * add the external statements to debug.h
 *
 * all the programs have to #include debug.cc.
 */

#undef DEBUG

#define SUBL(s) s ## _debug_level
#define SUBLEVEL(s) SUBL(s)

#define TOSTRING(s) #s
#define STR(s) TOSTRING(s)

#ifdef NDEBUG
#define DEBUG(level, s)
#else
#define DEBUG(level, s) \
  if (level <= ::debug_level + MODULE::debug_level + \
      MODULE::SUBLEVEL(SUBMODULE) ){ \
    std::cout << STR(MODULE) << ": " <<  s << std::endl; \
  }

// the global one
extern int debug_level;

namespace math
{
  extern int debug_level;
}

namespace interaction
{
  extern int debug_level;
  extern int interaction_debug_level;
  extern int pairlist_debug_level;
  extern int forcefield_debug_level;
}

namespace io
{
  extern int debug_level;
  extern int trajectory_debug_level;
  extern int input_debug_level;
  extern int topology_debug_level;
}

namespace simulation
{
  extern int debug_level;
  extern int simulation_debug_level;
  extern int topology_debug_level;
  extern int system_debug_level;
  extern int core_debug_level;
}

namespace algorithm
{
  extern int debug_level;
  extern int constraint_debug_level;
  extern int integration_debug_level;
  extern int algorithm_debug_level;
  extern int temperature_debug_level;
}
  
#endif
