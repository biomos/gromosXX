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
 * @file pscale.h
 * periodic scaling
 */

#ifndef INCLUDED_PSCALE_H
#define INCLUDED_PSCALE_H

namespace interaction
{
  struct dihedral_type_struct;
  
  /**
   * @class Periodic_Scaling
   * implements periodic scaling (of topology parameters)
   */
  class Periodic_Scaling : public interaction::Interaction
  {
  public:
    /**
     * Constructor.
     */
    Periodic_Scaling(interaction::Forcefield & ff,
		     simulation::Parameter const & param);
    
    /**
     * Destructor.
     */
    virtual ~Periodic_Scaling(){}
    
    /**
     * Periodic Scaling algorithm
     */
    virtual int calculate_interactions(topology::Topology &topo, 
				       configuration::Configuration &conf,
				       simulation::Simulation &sim);

    /**
     * initialises data structures.
     * adds additional types for the dihedral potentials,
     * sets up a map between j-value restraint and dihedral angle.
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    
  private:
    /**
     * scale a force constant
     */
    double scale(double t, double T, double s);


    /**
     * private dihedral types
     */
    std::vector<interaction::dihedral_type_struct>  m_dihedral_types;

  };
  
} // interaction

#endif

