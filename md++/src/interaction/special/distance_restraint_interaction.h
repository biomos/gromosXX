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
 * @file distance_restraint_interaction.h
 * distance restraining
 */

#ifndef INCLUDED_DISTANCE_RESTRAINT_INTERACTION_H
#define INCLUDED_DISTANCE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class distance_restraint_interaction
   * calculates the distance restraining interaction
   */
  class Distance_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Distance_Restraint_Interaction() : Interaction("DistanceRestraint"), 
        exponential_term(0.0) {}
    
    /**
     * Destructor.
     */
    virtual ~Distance_Restraint_Interaction() {}

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
    /**
     * exponential term used in averaging
     */
    double exponential_term;
    /**
     * Determine the dimensions for which the distance restraint applies
     *
     * Distance restraints may be applied to selected dimensions only by specifying the appropriate
     * code for RAH in the DISTANCERESSPEC or PERTDISRESSPEC block: 
     * 
     * @verbatim
                                       value of RAH
                      --------------------------------------------
Dimensions to         Half harmonic   Full harmonic  Half harmonic
apply restraint       repulsive                      attractive

 x, y, z              -1                0             1
 x, y                  9               10            11
 x, z                 19               20            21
 y, z                 29               30            31
 x                    39               40            41
 y                    49               50            51
 z                    59               60            61
    * @endverbatim
    */
    std::map<int,math::Vec> rah_map;

    
  };
  
} // interaction

#endif
