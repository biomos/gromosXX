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
 * @file energy_calculation.h
 * calculates the (total) energies and updates the averages
 */

#ifndef INCLUDED_ENERGY_CALCULATION_H
#define INCLUDED_ENERGY_CALCULATION_H

namespace algorithm
{
  /**
   * @class Energy_Calculation
   * calculates total energies, updates the averages
   */
  class Energy_Calculation : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Energy_Calculation() : Algorithm("EnergyCalculation") {}

    /**
     * Destructor.
     */
    virtual ~Energy_Calculation(){}
    
    /**
     * calculate the totals
     * update the averages
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

    // conf.old().energies.calculate_totals() (in apply()) reads
    // kinetic_energy, which Temperature_Calculation<gpuBackend> may
    // have only queued asynchronously on the GPU rather than finished
    // computing on the host yet -- see Algorithm::needs_finalized_gpu_
    // state()'s doc comment.
    virtual bool needs_finalized_gpu_state() const override { return true; }

    /**
     * apply() touches conf.old().energies (CPU-resident) and
     * conf.current().averages (plain CPU bookkeeping) -- confirmed by
     * reading the code, not guessing, same as Pressure_Calculation's
     * identical note. Narrowing away from the base class default
     * (MIRROR_ALL) still matters for the same reason documented
     * previously (avoids an unrelated full POS publish+invalidate every
     * step -- see PERFORMANCE.md's "Architecture direction" section),
     * but this class is no longer a true 0u: MIRROR_ENERGY is the one
     * real GPU-mirror-tracked field it needs. The GPU-native bonded/
     * special terms (CUDA_Angle_Interaction etc.) atomicAdd their
     * per-energy-group contributions into gpu::Configuration's energy_*
     * buffers instead of each keeping a private managed-memory scratch
     * buffer and touching it from the host every step (a real unified-
     * memory page-fault migration per touch -- see git history for the
     * profiling that found this); this is the single publish point that
     * copies them into conf.old().energies.{bond,angle,improper,
     * dihedral,posrest}_energy before calculate_totals() below reads
     * them. Same MIRROR_VIRIAL precedent as Pressure_Calculation.
     */
    virtual unsigned gpu_mirror_touches() const override { return gpu::MIRROR_ENERGY; }

  };
  
} // algorithm

#endif

