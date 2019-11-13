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
     * @param charge the charge
     */
    MM_Atom(unsigned index
          , double charge = 0.0
          , math::Vec pos = {0.0,0.0,0.0}
          , math::Vec force = {0.0,0.0,0.0}
          , double cos_charge = 0.0
          , math::Vec cosV = {0.0,0.0,0.0}
          ) :
              index(index)
            , charge(charge)
            , pos(pos)
            , force(force)
            , cos_charge(cos_charge)
            , cosV(cosV)
    {}
    
    /**
     * the index of this atom in the topology
     */
    unsigned int index;

    /**
     * the charge
     */
    mutable double charge;

    /**
     * the position
     */
    mutable math::Vec pos;

    /**
     * the force
     */
    mutable math::Vec force;

    /**
     * the charge-on-spring charge
     */
    mutable double cos_charge;

    /** 
     * charge-on-spring (distance vector between cos and real atom)
     */
    mutable math::Vec cosV;
    
    /**
     * less-than comparison operator
     */
    bool operator<(const MM_Atom & a) const {
      return index < a.index;
    }
  };
}
#endif	/* MM_ATOM_H */