/**
 * @file qm_atom.h
 * topological data structure for QM/MM
 */

#ifndef QM_ATOM_H
#define	QM_ATOM_H

namespace topology {
  /**
   * @struct qm_atom_struct
   * A structure to hold the atomic parameter for a quantum calculation
   */
  struct qm_atom_struct {
    /**
     * Constructor
     * @param index the index of the atom
     * @param atomic_number the atomic number of the atom, default hydrogen
     * @param link the atom is a link atom, default false
     * @param topology_charge the charge this atom has in the initial topology, default 0.0
     */
    qm_atom_struct(unsigned int index, unsigned int atomic_number = 1, bool link = false,
    double topology_charge = 0.0) :
    index(index), atomic_number(atomic_number), link(link),
    topology_charge(topology_charge) {}
/**
     * Copy constructor
     * @param a the atom to copy
     */
    qm_atom_struct(const qm_atom_struct & a) : index(a.index),
    atomic_number(a.atomic_number), link(a.link),
    topology_charge(a.topology_charge) {
    }
    /**
     * assignment operator
     * @param a the atom to assign
     * @return a reference to this atom
     */
    qm_atom_struct & operator=(const qm_atom_struct & a) {
      index = a.index; atomic_number = a.atomic_number; link = a.link;
      topology_charge = a.topology_charge;
      return *this;
    }
    /**
     * index pointer of this atom
     */
    unsigned int index;
    /**
     * atomic number of the atom i.e. 1 for H, 2 for He, etc...
     */
    unsigned int atomic_number;
    /**
     * atom is a link atom
     */
    bool link;
    /**
     * The charge this atom used to have in the initial topology
     */
    double topology_charge;
    /**
     * less-than comparison operator
     */
    bool operator<(const qm_atom_struct & rhs) const {
      return index < rhs.index;
    }
  };
}

#endif	/* QM_ATOM_H */

