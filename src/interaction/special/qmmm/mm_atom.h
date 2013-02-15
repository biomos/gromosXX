/**
 * @file mm_atom.h structure to store MM atom information
 */

#ifndef MM_ATOM_H
#define	MM_ATOM_H

namespace interaction {
/**
   * @struct MM_Atom
   * A structure to hold information of MM atoms, i.e. point charges
   */
  struct MM_Atom {
    /**
     * Constructor
     * @param index index of MM atom in topology
     * @param pos the position
     * @param charge the charge
     */
    MM_Atom(unsigned int index, const math::Vec & pos, double charge) :
    index(index), pos(pos), charge(charge) {
    }
    /**
     * the index of this atom in the topology
     */
    unsigned int index;
    /**
     * the position
     */
    math::Vec pos;
    /**
     * the charge
     */
    double charge;
  };
  
  /**
   * function to determine the MM atoms to include as point charges
   * based on a cutoff criterion
   * @param mm_atoms the MM atoms determined
   */
  void determine_mm_atoms(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation & sim,
            std::vector<MM_Atom> & mm_atoms);
}

#endif	/* MM_ATOM_H */

