/**
 * @file src/interaction/interaction.h
 * gathers common include directives for the interaction library.
 */

namespace interaction
{
  enum virial_enum { no_virial = 0, molecular_virial = 1, atomic_virial = 2 };
}

// general
#include "interaction/interaction.h"
#include "forcefield/forcefield.h"

// nonbonded base
#include "interaction/nonbonded/storage.h"
#include "interaction/nonbonded/nonbonded_base.h"

// nonbonded filter
#include "pairlist/filter/filter.h"
#include "pairlist/filter/exclusion_filter.h"
#include "pairlist/filter/perturbation_filter.h"
#include "pairlist/filter/range_filter.h"

// nonbonded pairlist algorithm
#include "pairlist/pairlist_algorithm.h"
#include "pairlist/standard_pairlist_algorithm.h"

// nonbonded pairlist
#include "pairlist/pairlist.h"

// nonbonded interaction
#include "interaction/nonbonded/nonbonded_innerloop.h"
#include "interaction/nonbonded/perturbed_nonbonded_innerloop.h"
#include "interaction/nonbonded/nonbonded_interaction.h"
#include "interaction/nonbonded/perturbed_nonbonded_interaction.h"

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

// and provide a general spec
#include "interaction/nonbonded/nonbonded_spec.h"

#ifndef NDEBUG
namespace interaction
{
  extern int debug_level;
  extern int interaction_debug_level;
  extern int pairlist_debug_level;
  extern int forcefield_debug_level;
}
#endif
