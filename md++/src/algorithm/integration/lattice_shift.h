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
 * @file lattice_shift.h
 * keeping track of lattice shifts
 */

#ifndef INCLUDED_LATTICE_SHIFT_H
#define INCLUDED_LATTICE_SHIFT_H

#include <type_traits>
namespace algorithm
{
  /**
   * @class Lattice_Shift_Tracker
   * keeps track of lattice shifts
   */
  template <typename Backend = util::cpuBackend>
  class Lattice_Shift_Tracker : public AlgorithmT<Backend>
  {
  public:
    /**
     * Constructor.
     */
    Lattice_Shift_Tracker() : AlgorithmT<Backend>("Lattice_Shift_Tracker") {
        static_assert(!std::is_same_v<Backend, util::gpuBackend>,
                      "Lattice_Shift_Tracker is not implemented for this backend.");
                  }

    /**
     * Destructor.
     */
    virtual ~Lattice_Shift_Tracker(){}
    
    /**
     * put CG into box and keep track of shift
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
    
  protected:
    template<math::boundary_enum b>
    void _apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);
  };
}
#endif	/* INCLUDED_LATTICE_SHIFT_H */

