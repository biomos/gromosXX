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
       * Constractor
       */
      Storage(){}

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
    std::vector<configuration::Energy> perturbed_energy_derivatives;
    /**
     * (longrange) virial storage.
     */
    math::Matrix virial_tensor;

    /**
     * zero all entities
     */
    void zero()
    {
      force = 0.0;
      energies.zero();

      for(size_t s=0, s_to = perturbed_energy_derivatives.size();
	  s != s_to; ++s)
	perturbed_energy_derivatives[s].zero();

      virial_tensor = 0.0;
    }

  };
  
} // interaction

#endif    
