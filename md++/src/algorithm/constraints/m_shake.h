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
 * @file m_shake.h
 * the m_shake algorithm.
 */

#ifndef INCLUDED_M_SHAKE_H
#define INCLUDED_M_SHAKE_H

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm
{
  /**
   * @class M_Shake
   * implements the m_shake algorithm.
   */
  class M_Shake : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    M_Shake(double const tolerance = 0.000001,
	  int const max_iterations = 1000,
	  std::string const name = "M_Shake");

    /**
     * Destructor.
     */
    virtual ~M_Shake();
        
    /**
     * apply m_shake.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

    /**
     * tolerance.
     */
    double const & tolerance()const {return m_tolerance;}
    /**
     * max iterations.
     */
    int const & max_iterations()const {return m_max_iterations;}

    /**
     * initialize startup positions and velocities
     * if required.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

  protected:

    /**
     * shake tolerance
     */
    double m_tolerance;
    /**
     * max iterations
     */
    const int m_max_iterations;
    /** 
     * rank and size for parallelization
     */
    int m_rank, m_size;
    /**
     * the factor matrix
     */
    math::Matrix factor;
    /**
     * the constraint lengths squared
     */
    math::Vec constr_length2;
    /**
     * inverted masses
     */
    math::Vec mass_i;

    inline
    int m_shake_molecule
    (
     configuration::Configuration & conf,
     bool & convergence,
     int first,
     double dt2i,
     const std::vector<topology::two_body_term_struct> & constr,
     math::GenericVec<math::Vec> const & dist_old,
     bool do_virial
     );
    
    void solvent
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     double dt, int const max_iterations,
     int & error
     );

  };
  
} //algorithm

#endif
