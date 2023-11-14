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
     * eds endstates (longrange) force storage.
     */
    std::vector<math::VArray> force_endstates;
    /**
     * multiAEDS endstates (longrange) force storage.
     */
    std::map<std::vector<int>, math::VArray> force_mult_endstates;
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
    /**
     * eds endstates (longrange) virial storage.
     */
    std::vector<math::Matrix> virial_tensor_endstates;
    /**
     * multi eds endstates (longrange) virial storage.
     */
    std::map<std::vector<int>, math::Matrix> virial_tensor_mult_endstates;
     /**
      * (longrange) electric field storage.
      */
    math::VArray electric_field;
    /**
     * indices of domain decomposition
     */
    std::vector<int> domain;
    /**
     * group wise forces
     */
    std::vector<std::vector<math::VArray> > force_groups;

    /**
     * zero all entities
     */
    void zero()
    {
      force = 0.0;
      electric_field = 0.0;
      energies.zero();
      perturbed_energy_derivatives.zero();
      virial_tensor = 0.0;
      assert(force_endstates.size() == virial_tensor_endstates.size());
      unsigned int size = force_endstates.size();
      for(unsigned int i = 0; i< size; i++){
        force_endstates[i] = 0.0;
        virial_tensor_endstates[i] = 0.0;
      }
      domain.clear();
      size = force_groups.size();
      for(unsigned int i = 0; i < size; ++i) {
        for(unsigned int j = 0; j < size; ++j) {
          force_groups[i][j] = 0.0;
        }
      }
      //MULTIAEDS
      assert(virial_tensor_mult_endstates.size() == force_mult_endstates.size());
      for(auto i: virial_tensor_mult_endstates){
	      virial_tensor_mult_endstates[i.first] = 0.0;
        force_mult_endstates[i.first] = 0.0;
      }
    }

  };
  
} // interaction

#endif    
