/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
     * ORIOL_GAMD
     * GAMD (longrange) force storage.
     */
    std::vector<math::VArray> force_gamd;
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
     * GAMD (longrange) virial storage.
     */
   std::vector<math::Matrix> virial_tensor_gamd;
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
      //ORIOL_GAMD
      assert(force_gamd.size() == virial_tensor_gamd.size());
      unsigned int gamdsize = force_gamd.size();
      for(unsigned int i = 0; i< gamdsize; i++){
           force_gamd[i] = 0.0;
           virial_tensor_gamd[i] = 0.0;
      }
      domain.clear();
      size = force_groups.size();
      for(unsigned int i = 0; i < size; ++i) {
        for(unsigned int j = 0; j < size; ++j) {
          force_groups[i][j] = 0.0;
        }
      }
    }

  };
  
} // interaction

#endif    
