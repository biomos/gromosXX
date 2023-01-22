/**
 * @file mm_atom.h
 * MM atom structure for QM/MM
 */

#ifndef INCLUDED_MM_ATOM_H
#define	INCLUDED_MM_ATOM_H

namespace interaction {
/**
   * @struct MM_Atom
   * A structure to hold information of MM atoms
   */
  struct MM_Atom {
    /**
     * Constructor
     * @param index index of MM atom in topology
     * @param pos the position
     * @param atomic_number - the atomic number of the atom
     * @param charge the charge
     */
    MM_Atom(unsigned index
          , math::Vec pos = {0.0,0.0,0.0}
          , unsigned atomic_number = 0
          , double charge = 0.0
          , math::Vec force = {0.0,0.0,0.0}
          , double cos_charge = 0.0
          , math::Vec cosV = {0.0,0.0,0.0}
          , math::Vec cos_force = {0.0,0.0,0.0}
          , bool polarisable = false
          ) :
              index(index)
            , atomic_number(atomic_number)
            , pos(pos)
            , charge(charge)
            , force(force)
            , cos_charge(cos_charge)
            , cosV(cosV)
            , cos_force(cos_force)
            , is_polarisable(polarisable)
    {}
    
    /**
     * the index of this atom in the topology
     */
    unsigned int index;

    /**
     * the atomic number of the atom
     */
    unsigned int atomic_number;

    /**
     * the position
     */
    mutable math::Vec pos;

    /**
     * the charge
     */
    mutable double charge;

    /**
     * the force
     */
    mutable math::Vec force;

    /**
     * the charge-on-spring charge
     */
    mutable double cos_charge;

    /** 
     * charge-on-spring vector (pointing from atom to cos)
     */
    mutable math::Vec cosV;

    /** 
     * charge-on-spring force
     */
    mutable math::Vec cos_force;

    /** 
     * charge-on-spring force
     */
    mutable bool is_polarisable;
    
    /**
     * less-than comparison operator
     */
    bool operator<(const MM_Atom & a) const {
      return index < a.index;
    }
  };
}
#endif	/* MM_ATOM_H */