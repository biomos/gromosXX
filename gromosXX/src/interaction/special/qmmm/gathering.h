/**
 * @file gathering.h gathering for QM/MM simulations
 */

#ifndef GATHERING_H
#define	GATHERING_H

namespace interaction {
  /**
   * function to gather the atoms in the QM zone.
   */
  void gather_qmzone(topology::Topology& topo,
          configuration::Configuration& conf,
          simulation::Simulation& sim,
          math::VArray & qm_pos);

  /**
   * function to gather the MM atoms with respect to
   * the atoms in the QM zone.
   */
  void gather_mm_atoms(topology::Topology& topo,
          configuration::Configuration& conf,
          simulation::Simulation& sim,
          math::VArray & qm_pos,
          std::vector<interaction::MM_Atom>& mm_atoms);
}


#endif	/* GATHERING_H */

