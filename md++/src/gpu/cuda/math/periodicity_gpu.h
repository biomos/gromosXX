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
 * @file periodicity_gpu.h
 * implementation of the periodic boundary condition functions for GPU.
 */

#pragma once

#undef MODULE
#undef SUBMODULE
#define MODULE gpu
#define SUBMODULE math

// #include "gpu/cuda/memory/topology_struct.h"
// #include "gpu/cuda/memory/configuration_struct.h"
// #include "math/periodicity.h"
// #include "gpu/cuda/kernels/hello_world.h"

namespace configuration {
  class Configuration;
}

namespace topology {
  class Topology;
}

namespace gpu
{
  /**
   * @class PeriodicityGpu
   * the periodic boundary condition functions for GPU.
   */
  template<math::boundary_enum b>
  class PeriodicityGpu : public math::Boundary_Implementation<b>
  {
  public:
    /**
     * Constructor.
     * If b is any no specific code will be generated,
     * otherwise one can specify optimized code.
     */
    PeriodicityGpu(math::Box const & bb);
    /**
     * puts a vector into the box (centered at (0, 0, 0).
     */
    template<typename VecType>
    void put_into_box(VecType &v)const;
    /**
     * puts a vector into the box centered at (Kx/2, Ly/2, Mz/2).
     */
    void put_into_positive_box(math::Vec &v)const;
    /**
     * put chargegroups into the box.
     */
    void put_chargegroups_into_box(configuration::Configuration & conf,
				   topology::Topology const & topo )const;
    
    /**
     * put chargegroups into the box and save the lattice shifts.
     */
    void put_chargegroups_into_box_saving_shifts(
                                   configuration::Configuration & conf,
				   topology::Topology const & topo )const;


    void gather_chargegroups(configuration::Configuration & conf, 
			     topology::Topology const & topo)const;

    void gather_molecules_into_box(configuration::Configuration & conf, 
				   topology::Topology const & topo)const;
    
  };
  
} // math