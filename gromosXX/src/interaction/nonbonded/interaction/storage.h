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
     * (longrange) force storage.
     */
    math::VArray force;
    /**
     * (longrange) energy storage.
     */
    configuration::Energy energies;
    /**
     * (longrange) potential energy lambda derivative storage.
     */
    configuration::Energy perturbed_energy_derivatives;
    /**
     * (longrange) virial storage.
     */
    math::Matrix virial_tensor;
    
  };
  
} // interaction

#endif    
