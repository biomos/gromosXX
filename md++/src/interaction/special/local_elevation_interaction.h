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
 * @file local_elevation_interaction.h
 * local elevation umbrella samping interaction
 */

#ifndef INCLUDED_LOCAL_ELEVATION_INTERACTION_H
#define	INCLUDED_LOCAL_ELEVATION_INTERACTION_H

namespace interaction {

  /**
   * @class xray_restraint_interaction
   * calculates the xray restraining interaction
   */ class Local_Elevation_Interaction : public Interaction {
  public:

    /**
     * Constructor.
     */
    Local_Elevation_Interaction() : Interaction("Local Elevation") {}
    /**
     * Destructor.
     */
    virtual ~Local_Elevation_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os = std::cout,
            bool quiet = false);
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);
  protected:
  };
} // interaction
#endif	/* INCLUDED_LOCAL_ELEVATION_INTERACTION_H */

