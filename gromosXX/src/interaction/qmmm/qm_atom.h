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
 * @file qm_atom.h
 * QM atom structure for QM/MM
 */

#ifndef INCLUDED_QM_ATOM_H
#define	INCLUDED_QM_ATOM_H

namespace interaction {
  /**
   * @struct QM_Atom
   * A structure holding information on QM atoms
   */
  struct QM_Atom {
    /**
     * Constructor
     * @param index - the index of the atom in topology
     * @param atomic_number - the atomic number of the atom
     * @param pos - 3D vector of atomic position
     * @param linked - the index of atom it is linked to - default is 0 (none)
     */
    QM_Atom(unsigned index
         , math::Vec pos = {0.0,0.0,0.0}
         , unsigned atomic_number = 1
         , bool is_linked = false
         ) : pos(pos)
           , force(0.0)
           , qm_charge(0.0)
           , index(index)
           , atomic_number(atomic_number)
           , is_linked(is_linked)
    {}
    
    /**
     * Copy constructor
     */
    QM_Atom(const QM_Atom & a) : pos(a.pos)
                               , force(a.force)
                               , qm_charge(a.qm_charge)
                               , index(a.index)
                               , atomic_number(a.atomic_number)
                               , is_linked(a.is_linked)
    {}

    /**
     * The coordinate
     */
    mutable math::Vec pos;

    /**
     * The force
     */
    mutable math::Vec force;

    /**
     * The charge calculated from QM
     */
    mutable double qm_charge;

    /**
     * the index of the atom in topology - starts with 0
     */
    const unsigned index;
    
    /**
     * the atomic number of the atom
     */
    const unsigned atomic_number;
    
    /**
     * if atom is linked to MM atom
     */
    mutable bool is_linked;

    /**
     * less-than comparison operator
     */
    bool operator<(const QM_Atom & a) const {
      return index < a.index;
    }
  };
}

#endif	/* QM_ATOM_H */

