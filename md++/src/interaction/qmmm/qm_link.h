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
 * @file qm_link.h
 * QM link class - implements QM-MM link atom schemes
 */

#ifndef INCLUDED_QM_LINK_H
#define	INCLUDED_QM_LINK_H

namespace interaction {
  struct QM_Atom;
  /**
   * @class QM_Link
   * Bonds between QM and MM atoms
   */
  struct QM_Link {
    /**
     * Constructor
     * @param cap_atom Capping atom (QM_Atom)
     * @param qm_index Index of QM atom
     * @param mm_index Index of MM atom
     */
    QM_Link(interaction::QM_Atom cap_atom
          , unsigned qm_index
          , unsigned mm_index
          ) : 
              cap_atom(cap_atom)
            , qm_index(qm_index)
            , mm_index(mm_index)
            , atomic_number(this->cap_atom.atomic_number)
            , pos(this->cap_atom.pos)
            , force(this->cap_atom.force)
            , qm_charge(this->cap_atom.qm_charge)
    {}

    /**
     * QM capping atom
     */
    mutable interaction::QM_Atom cap_atom;

    /**
     * QM link atom index
     */
    const unsigned qm_index;

    /**
     * MM link atom index
     */
    const unsigned mm_index;

    /**
     * Reference to capping atom atomic number (direct access)
     */
    const unsigned& atomic_number;

    /**
     * Reference to capping atom position (direct access)
     */
    math::Vec& pos;

    /**
     * Reference to capping atom force (direct access)
     */
    math::Vec& force;

    /**
     * Reference to capping atom charge (direct access)
     */
    double& qm_charge;

    /**
     * Update QM capping atom position
     */
    void update_cap_position(const math::Vec& qm_pos
                           , const math::Vec& mm_pos
                           , const double cap_length) const;

    /**
     * Distribute force on capping atom between QM and MM atom
     */
    void distribute_force(const math::Vec &qm_pos
                        , const math::Vec &mm_pos
                        , math::Vec &qm_force
                        , math::Vec &mm_force) const;

    /**
     * less-than comparison operator
     * Sorting in set by QM atom index, then by MM atom index
     * Allows multiple QM-MM links
     */
    bool operator<(const QM_Link& l) const {
      return qm_index < l.qm_index ||
            (qm_index == l.qm_index && mm_index < l.mm_index);
    }
  };
}

#endif	/* QM_LINK_H */

