/**
 * @file debug.h
 * define debug macros.
 * @page debug Debug Levels
 *  -# : no output
 *  -# : important program decisions
 *  -# : important program path
 *  -# : 
 *  -# : program path (not unimportant functions)
 *  -# : function outline
 *  -# :
 *  -# : function detail
 *  -# :
 *  -# :
 *  -# : everything
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
#define DEBUG
#else
#define DEBUG(level, s) \
  if (level <= ::debug_level + MODULE::debug_level + \
      MODULE::SUBLEVEL(SUBMODULE) ){ \
    std::cout << STR(MODULE) << ": " <<  s << std::endl; \
  }

// the global one
extern int debug_level;

namespace interaction
{
  extern int debug_level;
  extern int interaction_debug_level;
}

namespace io
{
  extern int debug_level;
  extern int trajectory_debug_level;
  extern int input_debug_level;
}

namespace simulation
{
  extern int debug_level;
  extern int simulation_debug_level;
  extern int topology_debug_level;
  extern int system_debug_level;
}

#endif
