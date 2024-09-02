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
 * @file print_block.h
 * routines to print out the various blocks.
 */

#ifndef INCLUDED_PRINT_BLOCK_H
#define INCLUDED_PRINT_BLOCK_H

namespace util
{
  struct Replica_Data;
}

namespace io
{
  /**
   * Print the dof coupling table from the MULTIBATH block.
   */
  void print_MULTIBATH_COUPLING(std::ostream &os, simulation::Multibath const &bath);

  /** 
   * Print the MULTIBATH block.
   */
  void print_MULTIBATH(std::ostream &os, simulation::Multibath const &bath,
		       configuration::Energy const &energy,
		       std::string title="TEMPERATURES");

  /** 
   * Print the DEGREES OF FREEDOM block.
   */
  void print_DEGREESOFFREEDOM(std::ostream &os, simulation::Multibath const &bath);

  /**
   * Print the PCOUPLE block.
   */
  void print_PCOUPLE(std::ostream &os, bool calc, 
		     math::pressure_scale_enum scale,
		     math::Matrix pres0, double comp, 
		     double tau, math::virial_enum vir,
                     int x_semi, int y_semi, int z_semi);

  /**
   * Print the PRESSURE block.
   */
  void print_PRESSURE(std::ostream &os, configuration::Configuration const & conf);
  
  /**
   * Print the ENERGY block.
   */
  void print_ENERGY(std::ostream &os, configuration::Energy const &e,
		    std::vector<unsigned int> const & energy_groups,
		    std::string const title = "ENERGIES",
		    std::string const type = "E_");
  
  /**
   * Print a matrix.
   */
  void print_MATRIX(std::ostream &os, math::Matrix const &m,
		    std::string const title);

  /**
   * Print the TIMESTEP block.
   */
  void print_TIMESTEP(std::ostream &os, double const steps, double const time);

  /**
   * Print the CENTREOFMASS block.
   */
  void print_CENTREOFMASS(std::ostream &os, double const ekin_trans, double const ekin_rot);
			    
   /** 
    * print RAMD information
    */
   void print_RAMD(std::ostream &os, configuration::Configuration const & conf, double lambda);

  /**
   * Print the SASA block.
   */
  void print_sasa(std::ostream &os, topology::Topology const & topo,
                    configuration::Configuration const & conf,
                    simulation::Simulation const & sim,
		    std::string const title);
  /**
   * Print sasa and volume averages and fluctuation.
   */
  void print_sasa_avg(std::ostream &os, topology::Topology const & topo,
                      std::vector<double> const &s,
                      std::vector<double> const &v, double &stot, double &vtot,
                      std::string const title, int volume);
  
  void print_sasa_fluct(std::ostream &os, topology::Topology const & topo,
                        std::vector<double> const &s,
                        std::vector<double> const &v, double &stot, double &vtot,
                        std::string const title, int volume);

} // io

#endif
