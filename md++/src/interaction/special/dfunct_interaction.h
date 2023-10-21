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
 * @file dfunct_interaction.h
 * dfunct
 */

#ifndef INCLUDED_DFUNCT_INTERACTION_H
#define INCLUDED_DFUNCT_INTERACTION_H

namespace interaction {

  /*
   * The DFUNCT interaction was described and used
   * to model transition states in this publication:
   * https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.0c01112
  */
  class DFunct_Interaction : public Interaction {

  public:

    /**
     * Constructor
     * 
     */
    DFunct_Interaction() : Interaction("DFunct") {}

    /**
     * Destructor
     * 
     */
    virtual ~DFunct_Interaction() = default;

    /**
     * Initializes the DFUNCT interaction
     * 
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param os Output stream
     * @param quiet Verbosity
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(topology::Topology& topo,
		                 configuration::Configuration& conf,
		                 simulation::Simulation& sim,
		                 std::ostream& os = std::cout,
		                 bool quiet = false);

    /**
     * Calculates the DFUNCT interaction
     * 
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @return 0 if successful, non-zero on failure
     */
    virtual int calculate_interactions(topology::Topology& topo,
				                               configuration::Configuration& conf,
				                               simulation::Simulation& sim);

  };

} // interaction

#endif /* INCLUDED_DFUNCT_INTERACTION_H */