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
 * @file lincs.h
 * the LINCS algorithm.
 */

#ifndef INCLUDED_LINCS_H
#define INCLUDED_LINCS_H

namespace algorithm {

  /**
   * Builds `lincs.coupled_constr`/`coef`/`sdiag` from a constraint
   * list (topo.solute().distance_constraints() or
   * topo.solvent(i).distance_constraints()) -- shared by
   * Lincs::init() and CUDA_Lincs::init() (cuda_lincs.cc), since the
   * setup itself is plain topology math with no GPU/CPU-specific
   * content.
   */
  void setup_lincs(topology::Topology const & topo,
                    topology::Compound::lincs_struct & lincs,
                    std::vector<topology::two_body_term_struct> const & constr,
                    unsigned int offset = 0);

  /**
   * @class Lincs
   * implements the lincs algorithm.
   */
  class Lincs : public Algorithm {
  public:
    /**
     * Constructor.
     */
    Lincs();

    /**
     * Destructor.
     */
    virtual ~Lincs();

    virtual int apply(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);


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

    /**
     * initialize startup positions and velocities
     * if required.
     */
    int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            std::ostream & os = std::cout,
            bool quiet = false);

  protected:

    /**
     * the atoms that are involved in the contraints
     */
    std::set<unsigned int> m_constrained_atoms;

  };

} //algorithm

#endif
