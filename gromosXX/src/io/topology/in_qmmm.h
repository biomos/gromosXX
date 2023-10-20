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
 * @file in_qmmm.h
 * read in a QM/MM specification file
 */
/**
 * @page qmmm QM/MM specification format
 * @date 22-01-2020
 *
 * A QM/MM specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref qmzone
 * - @ref qmunit
 */

#ifndef INCLUDED_IN_QMMM_H
#define INCLUDED_IN_QMMM_H

#include "../instream.h"

namespace io {

  /**
   * @class In_QMMM
   * reads in a QM/MM specification file
   */
  class In_QMMM : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_QMMM() {}
    /**
     * Constructor.
     */
    In_QMMM(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a QM/MM specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
    /**
     * Read in QM/MM units conversion factors
     */
    void read_units(const simulation::Simulation & sim
          , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param);
    /**
     * Read file paths for QM/MM trajectory
     */
    void read_trajectory_files(const simulation::Simulation & sim, simulation::Parameter::qmmm_struct::qm_param_struct* qm_param);
    /**
     * Read the map of atomic numbers to element names
     */
    void read_elements(const topology::Topology& topo
    , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param);
    /**
     * Read the map of IAC codes to atomic numbers
     */
    void read_iac_elements(topology::Topology& topo
    , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param);
    /**
     * Read the list of QM atoms
     */
    void read_zone(topology::Topology& topo
                    , simulation::Simulation& sim
                    , const std::string& blockname);
    /**
     * helper function to remove constraints from QM atoms
     */
    void remove_constraints(topology::Topology& topo);
  };
} // io

#endif
