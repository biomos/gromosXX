/**
 * @file src/simulation/core.h
 * gathers include directives for the core simulation headers.
 * these should be included first after math...
 */

#include "core/parameter.h"
#include "core/atom_iterator.h"
#include "core/atom_group_iterator.h"
#include "core/chargegroup_iterator.h"
#include "core/molecule_iterator.h"

#ifndef NDEBUG
namespace simulation
{
  extern int core_debug_level;
}

#endif
  
