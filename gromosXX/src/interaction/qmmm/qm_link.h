/**
 * @file qm_link.h
 * QM link class - implements QM-MM link atom schemes
 */

#ifndef INCLUDED_QM_LINK_H
#define	INCLUDED_QM_LINK_H

namespace interaction {
  class QM_Atom;
  class QM_Zone;
  /**
   * @class QM_Link
   * Bonds between QM and MM atoms
   */
  class QM_Link {
  public:
    /**
     * Constructor
     * @param Capping atom (QM_Atom)
     * @param Index of QM atom
     * @param Index of MM atom
     */
    QM_Link(unsigned qm_index
          , unsigned mm_index
          , interaction::QM_Atom cap_atom
          ) : 
              qm_index(qm_index)
            , mm_index(mm_index)
            , cap_atom(cap_atom)
            , atomic_number(this->cap_atom.atomic_number)
            , pos(this->cap_atom.pos)
            , force(this->cap_atom.force)
    {}

    /**
     * Index (capping atom has no index)
     */
    const int index = -1;

    /**
     * QM link atom index
     */
    const unsigned qm_index;

    /**
     * MM link atom index
     */
    const unsigned mm_index;

    /**
     * QM capping atom
     */
    mutable interaction::QM_Atom cap_atom;

    /**
     * Reference to capping atom index (direct access)
     */
    //const unsigned& index;

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
     * Update QM capping atom position
     */
    void update_cap_position(const math::Vec& qm_pos, const math::Vec& mm_pos, const double cap_length) const;

    /**
     * Distribute force on capping atom between QM and MM atom
     */
    void distribute_force(interaction::QM_Zone& qm_zone) const;

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

