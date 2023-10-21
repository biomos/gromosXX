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
 * @file read_input.h
 * input parameters
 */

#ifndef INCLUDED_READ_INPUT_H
#define INCLUDED_READ_INPUT_H

namespace util
{
  struct Replica_Data;
}

namespace io
{
  /**
   * read all input files necessary for the simulation.
   * calls
   * - read_parameter
   * - read_topology
   * - read_special
   * - read_configuration
   */
  int read_input
  (
   io::Argument const &args,
   topology::Topology &topo,
   configuration::Configuration &conf,
   simulation::Simulation & sim,
   algorithm::Algorithm_Sequence &md_seq,
   std::ostream & os = std::cout,
   bool quiet = false
   );

  /**
   * read the rest of the input files necessary for the repex simulation.
   * calls
   * - read_topology
   * - read_special
   * - read_configuration
   */
  int read_input_repex
  (
   io::Argument const &args,
   topology::Topology &topo,
   configuration::Configuration &conf,
   simulation::Simulation & sim,
   algorithm::Algorithm_Sequence &md_seq,
   int rank=0,
   int totalNumberOfThreads=-1,
   std::ostream & os = std::cout,
   bool quiet = false
   );

  /**
   * read in simulation parameter
   */
  int read_parameter
  (
   io::Argument const &args,
   simulation::Simulation & sim,
   std::ostream & os = std::cout,
   bool quiet = false
   );

  /**
   * read in topology, perturbation topology, MD sequence
   */
  int read_topology
  (
   io::Argument const &args,
   topology::Topology &topo,
   simulation::Simulation & sim,
   algorithm::Algorithm_Sequence &md_seq,
   std::ostream & os = std::cout,
   bool quiet = false
   );

  /**
   * read in a configuration
   */
  int read_configuration
  (
   io::Argument const &args,
   topology::Topology &topo,
   configuration::Configuration &conf,
   simulation::Simulation & sim,
   std::ostream & os = std::cout,
   bool quiet = false
   );
}

#endif
