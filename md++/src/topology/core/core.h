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
