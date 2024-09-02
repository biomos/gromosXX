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
 * @file gpu_shake.h
 * m_shake algroithm on the gpu
 */
#ifndef _GPU_SHAKE_H
#define	_GPU_SHAKE_H

#ifndef HAVE_LIBCUDART
#define gpu_status void
#endif
#include <algorithm/constraints/gpu_shake_thread.h>

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm
{
  /**
   * @class GPU_Shake
   * implements the m_shake algorithm.
   */
  class GPU_Shake : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    GPU_Shake(double const tolerance = 0.000001,
	  int const max_iterations = 1000,
	  std::string const name = "GPU_Shake");

    /**
     * Destructor.
     */
    virtual ~GPU_Shake();

    /**
     * apply gpu_shake.
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

    void solvent
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     double dt, int const max_iterations,
     int & error
     );

    /**
     * Pointer to all the GPU_Shake_Threads
     */
    std::vector<GPU_Shake_Thread *> m_shake_set;


  };

} //algorithm

#endif	/* _GPU_SHAKE_H */

