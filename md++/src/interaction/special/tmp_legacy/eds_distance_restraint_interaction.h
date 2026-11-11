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
 * @file eds_distance_restraint_interaction.h
 * eds perturbed distance restraining
 */

#ifndef INCLUDED_EDS_DISTANCE_RESTRAINT_INTERACTION_H
#define INCLUDED_EDS_DISTANCE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class eds_distance_restraint_interaction
   * calculates the eds perturbed distance restraining interaction
   */
  class Eds_Distance_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Eds_Distance_Restraint_Interaction() : Interaction("EdsDistanceRestraint"), 
        exponential_term(0.0) {}
    
    /**
     * Destructor.
     */
    virtual ~Eds_Distance_Restraint_Interaction() {}

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
    
  };
  
} // interaction

#endif
