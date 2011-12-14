/**
 * @file debug.h
 * define debug macros.
 */

/**
 * @page debug Debug
 *
 * @section Verbosity
 * To get debug messages a verbosity arguments has to be specified.
 * It is possible to set specific verbosity levels for each different
 * module and submodule.
 * The verbosity specified for a module is added to the general verbosity,
 * the one for a submodule to the verbosity of its parent module.
 * Negative verbosities are allowed.
 * 
 * <span style="font-size:larger"><b>
 * @verbatim  @verb 5 interaction:2 interaction:pairlist:-5 @endverbatim
 * </b></span><br>
 * This generally enables all debug messages with a level of less or equal to 5.
 * In the interaction module, debug messages with a level of less or equal to 7 are printed,
 * but for the submodule interaction:pairlist only the ones with a level of less or equal to 2.
 *
 * @section deblevels Debug Levels
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
 * @section debmodules Debug Modules
 * - math
 * - interaction
 *   - forcefield
 *   - interaction
 *   - pairlist
 *   - bonded
 *   - nonbonded
 *   - latticesum
 *   - special
 * - io
 *   - configuration
 *   - parameter
 *   - topology
 *   - forcefield
 * - configuration
 *   - configuration
 *   - energy
 * - topology
 * - algorithm
 *   - algorithm
 *   - constraints
 *   - temperature
 *   - pressure
 *   - integration
 * - simulation
 * - util
 *   - replica
 *   - leus
 *
 * @section Coding
 * in every file MODULE and SUBMODULE have to be defined
 * 
 * <span style="font-size:larger"><b>
 @verbatim #define MODULE module_name
 #define SUBMODULE submodule_name @endverbatim
 * </b></span><br>
 * and debug.h has to be included.
 *
 * To print out a debug message (if the appropriate verbosity
 * is selected), just use:<br>
 * <span style="font-size:larger"><b>
 @verbatim  DEBUG(level, "message var=" << var); @endverbatim
 * </b></span>
 * <br>
 *
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
  extern int math_debug_level;
}

namespace interaction
{
  extern int debug_level;
  extern int forcefield_debug_level;
  extern int interaction_debug_level;
  extern int pairlist_debug_level;
  extern int filter_debug_level;
  extern int bonded_debug_level;
  extern int nonbonded_debug_level;
  extern int latticesum_debug_level;
  extern int special_debug_level;
}

namespace io
{
  extern int debug_level;
  extern int configuration_debug_level;
  extern int parameter_debug_level;
  extern int topology_debug_level;
  extern int forcefield_debug_level;
}

namespace configuration
{
  extern int debug_level;
  extern int configuration_debug_level;
  extern int energy_debug_level;
}

namespace topology
{
  extern int debug_level;
  extern int topology_debug_level;
}

namespace algorithm
{
  extern int debug_level;
  extern int algorithm_debug_level;
  extern int constraints_debug_level;
  extern int temperature_debug_level;
  extern int pressure_debug_level;
  extern int integration_debug_level;
}

namespace simulation
{
  extern int debug_level;
  extern int simulation_debug_level;
}

namespace util
{
  extern int debug_level;
  extern int util_debug_level;
  extern int replica_debug_level;
  extern int leus_debug_level;
  extern int bs_leus_debug_level;
}

  
#endif
