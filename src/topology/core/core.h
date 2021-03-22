/**
 * @file core.h
 * basic concepts of a topology.
 * this subspace includes:
 * - @sa topology::Atom_Iterator
 * - @sa topology::Atomgroup_Iterator
 * - @sa topology::Chargegroup_Iterator
 * - @sa topology::Molecule_Iterator
 * - interaction term defintion: @sa body_term.h
 * - @sa Compound
 */


#include "atom_iterator.h"
#include "atomgroup_iterator.h"
#include "chargegroup_iterator.h"
#include "molecule_iterator.h"
#include "temperaturegroup_iterator.h"
#include "pressuregroup_iterator.h"
#include "../../util/virtual_atom.h"
#include "body_term.h"
#include "compound.h"
#include "qm_atom.h"
