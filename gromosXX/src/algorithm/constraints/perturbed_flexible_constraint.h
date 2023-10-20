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
 * @file perturbed_flexible_constraint.h
 * the perturbed flexible constraint algorithm.
 */

#ifndef INCLUDED_PERTURBED_FLEXIBLE_CONSTRAINT_H
#define INCLUDED_PERTURBED_FLEXIBLE_CONSTRAINT_H

namespace algorithm
{
  /**
   * @class Perturbed_Flexible_Constraint
   * implements the flexible constraint algorithm 
   * for perturbed distance constraints.
   */
  class Perturbed_Flexible_Constraint : public Flexible_Constraint
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Flexible_Constraint
    (
     double tolerance = 0.000001,
     int max_iterations = 1000,
     interaction::Forcefield *ff = NULL
     );
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Flexible_Constraint();
    
    /**
     * apply perturbed flexible constraints.
     * also applies unperturbed constraints
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    
    /**
     * initialize startup positions and velocities
     * if required.
     */
    int init(topology::Topology & topo,
	     configuration::Configuration & conf,
	     simulation::Simulation & sim,
	     std::ostream & os = std::cout,
	     bool quiet = false);

    void calc_distance
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim
     );

    void solute
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     int & error
     );

    void _store_lengths
    (
     configuration::Configuration & conf
     );

  protected:
  private:
    
    std::vector<double> m_perturbed_flex_len;

    template<math::boundary_enum B, math::virial_enum V>
    int _iteration
    (
     topology::Topology &topo,
     configuration::Configuration & conf,
     bool & convergence,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     double dt,
     math::Periodicity<B> const & periodicity
     );

    template<math::boundary_enum B>
    void _calc_distance
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim
     );

    template<math::boundary_enum B, math::virial_enum V>
    void _solute
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     int & error
     );

  };
  
} //algorithm

#endif
