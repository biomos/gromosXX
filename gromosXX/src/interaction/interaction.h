/**
 * @file src/interaction/interaction.h
 * gathers common include directives for the interaction library.
 */

// general
#include "interaction/interaction.h"
#include "forcefield/forcefield.h"

// nonbonded base
#include "interaction/storage.h"
#include "interaction/nonbonded_base.h"
#include "interaction/nonbonded_inner_loop.h"
#include "interaction/perturbed_nonbonded_inner_loop.h"
#include "interaction/nonbonded_inner_loop_virial.h"
#include "interaction/perturbed_nonbonded_inner_loop_virial.h"

// nonbonded filter
#include "pairlist/basic_filter.h"
#include "pairlist/perturbation_filter.h"
#include "pairlist/twinrange_filter.h"
#include "pairlist/twinrange_chargegroup_filter.h"

// nonbonded pairlist algorithm
#include "pairlist/basic_pairlist_algorithm.h"
#include "pairlist/chargegroup_range_pairlist_algorithm.h"

// nonbonded pairlist
#include "pairlist/basic_pairlist.h"

// nonbonded interaction
#include "interaction/nonbonded_interaction.h"
#include "interaction/perturbed_nonbonded_interaction.h"

// bonded interactions
#include "interaction/quartic_bond_interaction.h"
#include "interaction/perturbed_quartic_bond_interaction.h"
#include "interaction/harmonic_bond_interaction.h"
#include "interaction/perturbed_harmonic_bond_interaction.h"
#include "interaction/angle_interaction.h"
#include "interaction/perturbed_angle_interaction.h"
#include "interaction/improper_dihedral_interaction.h"
#include "interaction/perturbed_improper_dihedral_interaction.h"
#include "interaction/dihedral_interaction.h"
#include "interaction/perturbed_dihedral_interaction.h"

#ifndef NDEBUG
namespace interaction
{
  extern int debug_level;
  extern int interaction_debug_level;
  extern int pairlist_debug_level;
  extern int forcefield_debug_level;
}
#endif
