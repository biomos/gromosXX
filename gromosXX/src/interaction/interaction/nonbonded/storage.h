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
    /**
     * force accessor.
     */
    math::VArray & force() {return m_force;}
    /**
     * energy accessor.
     */
    simulation::Energy & energies() {return m_energy;}
    /**
     * lambda derivative of the potential energy accessor.
     */
    simulation::Energy & lambda_energies() {return m_lambda_energy;}
    /**
     * virial accessor.
     */
    math::Matrix & virial() {return m_virial;}
    
  protected:
    /**
     * (longrange) force storage.
     */
    math::VArray m_force;
    /**
     * (longrange) energy storage.
     */
    simulation::Energy m_energy;
    /**
     * (longrange) potential energy lambda derivative storage.
     */
    simulation::Energy m_lambda_energy;
    /**
     * (longrange) virial storage.
     */
    math::Matrix m_virial;
    
  };
  
} // interaction

#endif    
