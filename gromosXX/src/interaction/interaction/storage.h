/**
 * @file storage.h
 * store forces, energies and virial
 * for the longrange
 * nonbonded interactions.
 */

#ifndef INCLUDED_STORAGE_H
#define INCLUDED_STORAGE_H

namespace interaction
{
  /**
   * @class Storage
   * stores longrange nonbonded properties.
   */
  class Storage
  {
  public:
    math::VArray & force() {return m_force;}
    simulation::Energy & energies() {return m_energy;}
    math::Matrix & virial() {return m_virial;}
    
  protected:
    math::VArray m_force;
    simulation::Energy m_energy;
    math::Matrix m_virial;
    
  };
  
} // interaction

#endif    
