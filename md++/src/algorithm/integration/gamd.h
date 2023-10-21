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

/*
 * File:   gamd.h
 * Author: Oriol
 **/

#ifndef GAMD_H
#define	GMAD_H
namespace algorithm
{
  /**
   * @class GAMD
   * implements GAMD.
   */
  class GAMD : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    GAMD() : Algorithm("GAMD"){}

    /**
     * Destructor.
     */
    virtual ~GAMD(){}
    
    /**
     * 
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

  /**
   * calculate the mean, std and Vmax Vmin
   */
  void calc_gamd_std_mean(double V, int step, double *Vmax, double *Vmin, double *Vmean, long double *M2, long double *sigmaV);
  /**
   * calculate E threshold and k0
   */
  int calc_gamd_E_K(simulation::gamd_thresh_enum Eform, double sigma0, double Vmax, double Vmin, double Vmean, long double sigmaV, double *k0, double *k, double *E);
  /**
   * calculate interaction factor between acceleration regions
   */
  void calc_interaction_factor(int accelerationgroup, int accelerationgroup2, double *interaction_factor);
  };
   
} // algorithm

#endif
