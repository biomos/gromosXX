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
 * @file gpu_settle.h
 * Similar to settle, but running on the GPU
 */

#ifndef _GPU_SETTLE_H
#define	_GPU_SETTLE_H

// #ifdef HAVE_LIBCUDART
// #include "cukernel/cudaKernel.h"
// #endif

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm {
  class GPU_Settle : public Algorithm {
  public:
        /**
     * Constructor.
     */
    GPU_Settle() : Algorithm("GPU_Settle") {
    }

    /**
     * Destructor.
     */
    virtual ~GPU_Settle() {
    }

    /**
     * apply the SETTLE algorithm
     */
    virtual int apply(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);

    /**
     * initialize startup positions and velocities
     * if required.
     */
    virtual int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            std::ostream & os = std::cout,
            bool quiet = false);

    /**
     * accessor to the constrained atoms
     */
    std::set<unsigned int> & constrained_atoms() {
      return m_constrained_atoms;
    }

    /**
     * accessor to the constrained atoms
     */
    const std::set<unsigned int> & constrained_atoms() const {
      return m_constrained_atoms;
    }

  protected:
    /**
     * print an error message to std::cout if GPU_SETTLE fails
     */
    void printError(topology::Topology const & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            unsigned int mol);

    /** 
     * apply it for the solvent with correct periodic boundary conditions
     */
    void solvent(topology::Topology const & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            int & error);

    /**
     * the atoms that are involved in the contraints
     */
    std::set<unsigned int> m_constrained_atoms;
    /** 
     * rank and size for parallelization
     */
    int m_rank, m_size;
    /**
     * informations for the GPU
     */
    //gpu_status * gpu_stat;
  };
  
} // algorithm


#endif	/* _GPU_SETTLE_H */

