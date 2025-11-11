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
#include "AtomicRadii.h"

#include <cassert>
#include <cstdio>
#include <cmath>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/LJType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/GromosForceField.h"

using namespace gcore;

void utils::compute_atomic_radii_vdw(gcore::System & sys, const gcore::GromosForceField & gff) {
  static const double small = 1.0E-20;
  for (int m = 0; m < sys.numMolecules(); ++m) {
    for (int a = 0; a < sys.mol(m).topology().numAtoms(); ++a) {
      int atom_iac = sys.mol(m).topology().atom(a).iac();
      LJType lj(gff.ljType(AtomPair(atom_iac, atom_iac)));
      if (lj.c6() >= small) {
        sys.mol(m).topology().atom(a).setradius(0.5*(exp(log((2.0 * lj.c12()) / lj.c6()) / 6.0)));
      } else {
        sys.mol(m).topology().atom(a).setradius(0.0);
      }
    }
  }
  for(int s = 0; s < sys.numSolvents(); ++s) {
    for(int a = 0; a < sys.sol(s).topology().numAtoms(); ++a) {
      int atom_iac = sys.sol(s).topology().atom(a).iac();
      LJType lj(gff.ljType(AtomPair(atom_iac, atom_iac)));
      if (lj.c6() >= small) {
        sys.sol(s).topology().atom(a).setradius(0.5*(exp(log((2.0 * lj.c12()) / lj.c6()) / 6.0)));
      } else {
        sys.sol(s).topology().atom(a).setradius(0.0);
      }
    }
  }
}

void utils::compute_atomic_radii_vdw(int probe_iac, double probe_radius, gcore::System & sys, const gcore::GromosForceField & gff) {
  static const double small = 1.0E-20;
  for (int m = 0; m < sys.numMolecules(); ++m) {
    for (int a = 0; a < sys.mol(m).topology().numAtoms(); ++a) {
      int atom_iac = sys.mol(m).topology().atom(a).iac();
      LJType lj(gff.ljType(AtomPair(atom_iac, probe_iac)));
      if (lj.c6() >= small) {
        sys.mol(m).topology().atom(a).setradius((exp(log(2.0 * lj.c12() / lj.c6()) / 6.0)) - probe_radius);
      } else {
        sys.mol(m).topology().atom(a).setradius(0.0);
      }
    }
  }
  for(int s = 0; s < sys.numSolvents(); ++s) {
    for(int a = 0; a < sys.sol(s).topology().numAtoms(); ++a) {
      int atom_iac = sys.sol(s).topology().atom(a).iac();
      LJType lj(gff.ljType(AtomPair(atom_iac, probe_iac)));
      if (lj.c6() >= small) {
        sys.sol(s).topology().atom(a).setradius((exp(log(2.0 * lj.c12() / lj.c6()) / 6.0)) - probe_radius);
      } else {
        sys.sol(s).topology().atom(a).setradius(0.0);
      }
    }
  }
}
