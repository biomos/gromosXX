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
