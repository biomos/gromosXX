/**
 * @file src/interaction/interaction.h
 * gathers common include directives for the interaction library.
 */

#include "forcefield/parameter.h"
#include "interaction/interaction.h"
#include "forcefield/forcefield.h"

#include "pairlist/simple_pairlist.h"
#include "pairlist/twin_range_pairlist.h"
#include "pairlist/twin_range_pairlist_cg.h"

#include "interaction/nonbonded_interaction.h"
#include "interaction/nonbonded_virial_interaction.h"
#include "interaction/quartic_bond_interaction.h"
#include "interaction/harmonic_bond_interaction.h"
#include "interaction/angle_interaction.h"
#include "interaction/improper_dihedral_interaction.h"
#include "interaction/dihedral_interaction.h"

#ifndef NDEBUG
namespace interaction
{
  extern int debug_level;
  extern int interaction_debug_level;
  extern int pairlist_debug_level;
  extern int forcefield_debug_level;
}
#endif
