/**
 * @file io.h
 * gathers common include directives for io
 */

#include "argument.h"
#include "blockinput.h"
#include "message.h"
#include "GInStream.h"
#include "topology/InTopology.h"
#include "topology/InPerturbationTopology.h"
#include "trajectory/InTrajectory.h"
#include "trajectory/OutTrajectory.h"
#include "trajectory/OutG96Trajectory.h"
#include "trajectory/InFlexibleConstraints.h"
#include "input/InInput.h"
#include "print_block.h"

#ifndef NDEBUG

/**
 * @namespace io
 * provide the input/output routines.
 */
namespace io
{
  /**
   * the module debug level.
   */
  extern int debug_level;
  /**
   * debug level for the submodule trajectory.
   */
  extern int trajectory_debug_level;
  /**
   * debug level for the submodule input.
   */
  extern int input_debug_level;
  /**
   * debug level for the submodule topology.
   */
  extern int topology_debug_level;
}

#endif
