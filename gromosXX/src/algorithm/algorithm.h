/**
 * @file algorithm.h
 * gathers common include directives for algorithm
 */

#include "temperature/berendsen.h"
#include "pressure/berendsen.h"
#include "integration/leap_frog.h"
#include "constraint/shake.h"
#include "constraint/perturbed_shake.h"
#include "constraint/flexible_constraint.h"
#include "constraint/perturbed_flexible_constraint.h"

#include "algorithm/md_spec.h"
#include "algorithm/interaction_spec.h"

#include "algorithm/forcefield.h"
#include "algorithm/md_base.h"
#include "algorithm/md.h"
#include "algorithm/perturbation_md.h"


#ifndef NDEBUG
namespace algorithm
{
  extern int debug_level;
  extern int constraint_debug_level;
  extern int integration_debug_level;
  extern int algorithm_debug_level;
  extern int temperature_debug_level;
}
#endif
