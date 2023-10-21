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
 * @file ifp.h
 * interaction function parameter read-in interface
 */

#ifndef INCLUDED_IFP_H
#define INCLUDED_IFP_H

namespace io {

    /**
     * @class IFP
     * interface to read in interaction function parameters.
     */
    class IFP {
    public:

        /**
         * destructor
         */
        virtual ~IFP() {
        }

        /**
         * Read in the nonbonded interaction types (lennard-jones).
         */
        virtual void read_lj_parameter(std::vector<std::vector
                <interaction::lj_parameter_struct> >
                & lj_parameter,
                std::ostream & os = std::cout) = 0;

        /**
         * Read in the nonbonded interaction types (Coarse - grained lennard-jones).
         */
        virtual void read_cg_parameter(std::vector<std::vector
                <interaction::lj_parameter_struct> >
                & cg_parameter,
                std::ostream & os = std::cout) = 0;

        /**
         * Read in the nonbonded interaction types (SASA).
         */
        virtual void read_sasa_parameter(topology::Topology & topo,
                std::vector<topology::sasa_parameter_struct>
                & sasa_parameter) = 0;

    };

} // namespace io

#endif
