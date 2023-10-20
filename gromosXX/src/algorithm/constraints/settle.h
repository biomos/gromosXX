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
 * @file settle.h
 * the settle algorithm.
 */

#ifndef INCLUDED_SETTLE_H
#define INCLUDED_SETTLE_H

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm {

  /**
   * @class Settle
   * implements the settle algorithm for 3 site water models.
   *
   * Coordinates and velocities are reset according to
   * S. Miyamoto and P. A. Kollman, SETTLE: An Analytical Version of the SHAKE
   * and RATTLE Algorithm for Rigid Water Models, J. Comput. Chem 13, Issue 8, 
   * 1992, pp. 952-962
   *
   * The source code is annotaed according to the variables names used in the 
   * paper. 
   *
   * The constraint force and virial is not calculated as discribed in the paper
   * but from the displacement due to the constraint.
   */
  class Settle : public Algorithm {
  public:

    /**
     * Constructor.
     */
    Settle() : Algorithm("Settle") {
    }

    /**
     * Destructor.
     */
    virtual ~Settle() {
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
     * print an error message to std::cout if SETTLE fails
     */
    void printError(topology::Topology const & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            unsigned int atom, std::string message);

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
  };

} //algorithm

#endif
