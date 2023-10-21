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
 * @file perturbed_shake.h
 * the perturbed shake algorithm.
 */

#ifndef INCLUDED_PERTURBED_SHAKE_H
#define INCLUDED_PERTURBED_SHAKE_H

namespace algorithm
{
  
  /**
   * @class Perturbed_Shake
   * implements the shake algorithm for perturbed distance constraints.
   */
  class Perturbed_Shake : public Shake
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Shake(double const solute_tolerance = 0.000001,
        double const solvent_tolerance = 0.000001,
		    int const max_iterations = 1000);
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Shake();
        
    /**
     * apply perturbed SHAKE
     * (also calls SHAKE for the unperturbed bonds)
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

  protected:

    /**
     * do one iteration
     */
    template<math::boundary_enum B, math::virial_enum V>
    int perturbed_shake_iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     double tolerance,
     bool & convergence,
     int first,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     std::vector<topology::perturbed_two_body_term_struct>
     const & constr,
     double dt,
     math::Periodicity<B> const & periodicity,
     simulation::Simulation & sim); //ANITA

    /**
     * do a perturbed angle constraint iteration
     */
    template<math::boundary_enum B, math::virial_enum V>
    int perturbed_ang_constr_iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim,
     bool & convergence,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     math::Periodicity<B> const & periodicity
     );

    /**
     * do a perturbed dihedral constraint iteration
     */
    template<math::boundary_enum B, math::virial_enum V>
    int perturbed_dih_constr_iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim,
     bool & convergence,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     math::Periodicity<B> const & periodicity
     );

    /**
     * shake (perturbed) solute
     */
    template<math::boundary_enum B, math::virial_enum V>
    void perturbed_solute(topology::Topology const & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  int max_iterations,
			  int & error);
  // ANITA		  simulation::Simulation const & sim,
    
    // overwrites the other one, as g++ seems unable to compile this...!!!
    // seems to work with gcc 4.1 (by Clara & Nathan)
    /**
     * shake solvent (not perturbed)
     */
    /*
    template<math::boundary_enum B, math::virial_enum V>
    void solvent(topology::Topology const & topo,
		 configuration::Configuration & conf,
		 double dt, int const max_iterations,
		 int & error);*/
    
  };
  
} //algorithm

#endif
