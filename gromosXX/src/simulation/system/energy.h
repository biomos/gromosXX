/**
 * @file energy.h
 * storage of the energies.
 */

#ifndef INCLUDED_ENERGY_H
#define INCLUDED_ENERGY_H

namespace simulation
{
  /**
   * @class Energy
   * storage of energies.
   */
  class Energy
  {
  public:
    std::vector<double> bond_energy;
    std::vector<double> angle_energy;
    std::vector<double> improper_energy;
    std::vector<double> dihedral_energy;

    std::vector<std::vector<double> > lj_energy;
    std::vector<std::vector<double> > crf_energy;
    
    void zero();
    void resize(size_t s);
    
  };

} // simulation

#include "energy.tcc"

#endif
