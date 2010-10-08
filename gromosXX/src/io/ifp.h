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
         * Read in the harmonic bond parameter.
         */
        virtual void read_harmonic_bonds(std::vector<interaction::bond_type_struct> &b,
                std::ostream & os = std::cout) = 0;
        /**
         * Read in the quartic bond parameter.
         */
        virtual void read_g96_bonds(std::vector<interaction::bond_type_struct> &b,
                std::ostream & os = std::cout) = 0;
        /**
         * Read in the bond angle parameter.
         */
        virtual void read_angles(std::vector<interaction::angle_type_struct> &a,
                std::ostream & os = std::cout) = 0;

        /**
         * Read in the harmonic bond angle parameter.
         */
        virtual void read_harm_angles(std::vector<interaction::angle_type_struct> &a,
                std::ostream & os = std::cout) = 0;

        /**
         * Read in the improper dihedral parameter.
         */
        virtual void read_improper_dihedrals(std::vector<interaction::improper_dihedral_type_struct> &i,
                std::ostream & os = std::cout) = 0;

        /**
         * Read in the dihedral parameter.
         */
        virtual void read_dihedrals(std::vector<interaction::dihedral_type_struct> &d,
                std::ostream & os = std::cout) = 0;

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
