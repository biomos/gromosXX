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
         ) : index(index)
           , pos(pos)
           , atomic_number(atomic_number)
           , force(0.0)
           , qm_charge(0.0)
    {}
    /**
     * Copy constructor
     */
    QM_Atom(const QM_Atom & a) : index(a.index)
                               , pos(a.pos)
                               , atomic_number(a.atomic_number)
                               , force(a.force)
                               , qm_charge(a.qm_charge)
    {}

    /**
     * the index of the atom in topology - starts with 0
     */
    const unsigned index;

    /**
     * The coordinate
     */
    mutable math::Vec pos;
    
    /**
     * the atomic number of the atom
     */
    const unsigned atomic_number;

    /**
     * The force
     */
    mutable math::Vec force;

    /**
     * The charge calculated from QM
     */
    mutable double qm_charge;

    /**
     * less-than comparison operator
     */
    bool operator<(const QM_Atom & a) const {
      return index < a.index;
    }
  };
}

#endif	/* QM_ATOM_H */

