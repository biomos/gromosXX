/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
 *   - qmmm
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
 *   - timing
 *   - replica_exchange
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
#undef MPI_DEBUG

#define SUBL(s) s ## _debug_level
#define SUBLEVEL(s) SUBL(s)

#define TOSTRING(s) #s
#define STR(s) TOSTRING(s)

#ifdef NDEBUG
    #define DEBUG(level, s) ;
    #define MPI_DEBUG(level, s) ;

#else
    #define DEBUG(level, s) \
      if (level <= ::debug_level + MODULE::debug_level + \
          MODULE::SUBLEVEL(SUBMODULE) ){ \
        std::cout << STR(MODULE) << ": " <<  s << std::endl; \
      };

    #ifdef XXMPI
        #define MPI_DEBUG(level, s) \
          if (level <= ::debug_level + MODULE::debug_level + \
              MODULE::SUBLEVEL(SUBMODULE) ){ \
            MPI_Barrier(MPI_COMM_WORLD); \
            std::cout << STR(MODULE) << ": " <<  s << std::endl; \
            std::cout.flush(); \
            MPI_Barrier(MPI_COMM_WORLD);\
          };

    #else
        #define MPI_DEBUG(level, s)\
         DEBUG(level, s);

    #endif

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
      extern int qmmm_debug_level;
      extern int cuda_debug_level;
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
      extern int virtualatoms_debug_level;
    }

    namespace simulation
    {
      extern int debug_level;
      extern int simulation_debug_level;
    }

    namespace gpu
    {
      extern int debug_level;
      extern int kernel_debug_level;
      extern int constraints_debug_level;
      extern int interaction_debug_level;
      extern int pairlist_debug_level;
      extern int utils_debug_level;
    }

    namespace re
    {
        extern int debug_level;
        extern int replica_exchanger_debug_level;
        extern int replica_debug_level;
    }

    namespace util
    {
      extern int debug_level;
      extern int util_debug_level;
      extern int timing_debug_level;
      extern int leus_debug_level;
      extern int bs_leus_debug_level;
    }

#endif
